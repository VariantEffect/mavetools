from unittest import TestCase
import pandas as pd
from io import StringIO
from mavetools.mavedf.df_to_pandas import df_to_pandas


class Test(TestCase):

    # scenario 1 - valid input file
    def test_df_to_pandas(self):
        # typical case - return metadata
        data = StringIO(
            "# Accession:\n"
            "# Downloaded (UTC):\n"
            "# Licence:\n"
            "# Licence URL:\n"
            "accession,hgvs_nt,hgvs_pro,score\n"
            "urn:mavedb:00000011-a-1#27,c.[4=;5=;6=],p.Phe2=,0.0"
        )
        results = df_to_pandas(data, True)
        self.assertTrue(results[0].iat[0, 2] == "p.Phe2=")

        # typical case - do not return metadata
        data = StringIO(
            "# Accession:\n"
            "# Downloaded (UTC):\n"
            "# Licence:\n"
            "# Licence URL:\n"
            "accession,hgvs_nt,hgvs_pro,score\n"
            "urn:mavedb:00000011-a-1#27,c.[4=;5=;6=],p.Phe2=,0.0"
        )
        results = df_to_pandas(data, False)
        self.assertTrue(results.iat[0, 2] == "p.Phe2=")

        # edge case - ret_meta is incorrect type
        data = StringIO(
            "# Accession:\n"
            "# Downloaded (UTC):\n"
            "# Licence:\n"
            "# Licence URL:\n"
            "accession,hgvs_nt,hgvs_pro,score\n"
            "urn:mavedb:00000011-a-1#27,c.[4=;5=;6=],p.Phe2=,0.0"
        )
        with self.assertRaises(TypeError):
            df_to_pandas(data, "notabool")

    # scenario 2 - invalid filename
    def test_df_to_pandas_invalid_first_arg(self):
        # edge case - invalid filename
        with self.assertRaises(ValueError):
            df_to_pandas("notarealfile")

    # scenario 3 - missing arguments
    def test_df_to_pandas_missing_arguments(self):
        # edge case - missing first argument
        with self.assertRaises(TypeError):
            df_to_pandas()
