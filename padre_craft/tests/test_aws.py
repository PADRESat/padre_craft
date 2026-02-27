"""Tests util.py that interact with timestream"""

import os

import boto3
import pytest
from astropy import units as u
from astropy.timeseries import TimeSeries
from moto import mock_aws
from moto.core import DEFAULT_ACCOUNT_ID as ACCOUNT_ID
from moto.timestreamwrite.models import timestreamwrite_backends
from padre_meddea.housekeeping.calibration import get_calibration_func

import padre_craft.io.aws_db as aws
from padre_craft import _test_files_directory
from padre_craft.dirlist.dirlist import read_dirlist, summarize_dirlist
from padre_craft.io import read_file


@pytest.fixture(scope="function")
def aws_credentials():
    """Mocked AWS Credentials for moto."""
    os.environ["AWS_ACCESS_KEY_ID"] = "testing"
    os.environ["AWS_SECRET_ACCESS_KEY"] = "testing"
    os.environ["AWS_SECURITY_TOKEN"] = "testing"
    os.environ["AWS_SESSION_TOKEN"] = "testing"
    os.environ["AWS_DEFAULT_REGION"] = "us-east-1"


@pytest.fixture(scope="function")
def mocked_timestream(aws_credentials):
    """
    Return a mocked S3 client
    """
    with mock_aws():
        """Fixture to mock Timestream database and table."""
        client = boto3.client("timestream-write", region_name="us-east-1")
        database_name = "dev-padre_sdc_aws_logs"
        table_name = "dev-padre_measures_table"
        client.create_database(DatabaseName=database_name)

        client.create_table(
            DatabaseName=database_name,
            TableName=table_name,
            RetentionProperties={
                "MemoryStoreRetentionPeriodInHours": 24,
                "MagneticStoreRetentionPeriodInDays": 7,
            },
        )
        yield client


def test_record_timeseries_quantity_1col(mocked_timestream):
    timeseries_name = "meddea"
    ts = TimeSeries(
        time_start="2016-03-22T12:30:31",
        time_delta=3 * u.s,
        n_samples=5,
        meta={"name": timeseries_name},
    )
    print(os.getenv("SWXSOC_MISSION"))
    ts["FPTemp"] = [47984, 47994, 47884, 47984, 47984]

    aws.record_housekeeping(ts, "meddea")

    f = get_calibration_func("fp_temp")
    database_name = "dev-padre_sdc_aws_logs"
    table_name = "dev-padre_measures_table"

    backend = timestreamwrite_backends[ACCOUNT_ID]["us-east-1"]
    records = backend.databases[database_name].tables[table_name].records
    # Assert that there should be 5 records, one for each timestamp
    assert len(records) == len(ts["FPTemp"])

    for i, record in enumerate(records):
        # Assert the time is correct
        time = str(int(ts.time[i].to_datetime().timestamp() * 1000))
        assert record["Time"] == time
        assert record["MeasureName"] == timeseries_name
        # Check the MeasureValues
        measure_values = record["MeasureValues"]
        assert len(measure_values) == 1  # Only one column of data

        # Assert the measure name, value, and type
        temp4_measure = next(
            (mv for mv in measure_values if mv["Name"] == "fp_temp_deg_C"), None
        )
        assert temp4_measure is not None, "fp_temp_deg_C not found in MeasureValues"
        assert temp4_measure["Value"] == str(f(ts["FPTemp"][i]).value), (
            "MeasureValue does not match"
        )
        assert temp4_measure["Type"] == "DOUBLE", "MeasureValueType does not match"


def test_record_timeseries(mocked_timestream):
    hk_ts = read_file(
        _test_files_directory
        / "padre_get_MEDDEA_HOUSE_KEEPING_Data_1762493454480_1762611866270.csv"
    )
    aws.record_housekeeping(hk_ts, "meddea")
    database_name = "dev-padre_sdc_aws_logs"
    table_name = "dev-padre_measures_table"

    backend = timestreamwrite_backends[ACCOUNT_ID]["us-east-1"]
    records = backend.databases[database_name].tables[table_name].records
    # Assert that there should be 5 records, one for each timestamp
    assert (
        len(records) == len(hk_ts.time) - 3
    )  # removes 2 rows with zeros and one with future time


def test_record_dirlist(mocked_timestream):
    """Test recording dirlist summary to AWS Timestream."""
    # Read dirlist file and create summary
    test_dirlist_file = _test_files_directory / "padre_craft_dirlist_022426.txt"
    file_list = read_dirlist(test_dirlist_file)
    dirlist_summary = summarize_dirlist(file_list)

    # Record to AWS
    aws.record_dirlist(dirlist_summary)

    database_name = "dev-padre_sdc_aws_logs"
    table_name = "dev-padre_measures_table"

    backend = timestreamwrite_backends[ACCOUNT_ID]["us-east-1"]
    records = backend.databases[database_name].tables[table_name].records

    # Assert that there should be 1 record (one row in summary table)
    assert len(records) == len(dirlist_summary)
    assert len(records) == 1

    # Check the record details
    record = records[0]
    assert record["MeasureName"] == "dirlist"

    # Check that measure values contain expected columns
    measure_values = record["MeasureValues"]
    measure_names = [mv["Name"] for mv in measure_values]

    # Verify all expected columns are present
    expected_columns = [
        "file_count_meddea",
        "file_count_meddea_photon",
        "file_count_meddea_spectrum",
        "file_count_meddea_hk",
        "file_count_sharp",
        "file_count_total",
        "size_meddea",
        "size_meddea_photon",
        "size_meddea_spectrum",
        "size_meddea_hk",
        "size_sharp",
        "size_total",
    ]

    for col in expected_columns:
        assert col in measure_names, f"{col} not found in measure values"

    # Verify some specific values
    file_count_total = next(
        (mv for mv in measure_values if mv["Name"] == "file_count_total"), None
    )
    assert file_count_total is not None
    assert file_count_total["Value"] == str(dirlist_summary["file_count_total"][0])

    file_count_meddea = next(
        (mv for mv in measure_values if mv["Name"] == "file_count_meddea"), None
    )
    assert file_count_meddea is not None
    assert file_count_meddea["Value"] == str(dirlist_summary["file_count_meddea"][0])
