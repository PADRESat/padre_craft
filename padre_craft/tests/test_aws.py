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
