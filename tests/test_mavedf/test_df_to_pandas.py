from unittest import TestCase
import pandas as pd
from io import StringIO
from mavetools.mavedf.df_to_pandas import df_to_pandas


class Test(TestCase):

    # scenario 1 - call function on urn mavedb 00000001-a-1_scores

    # typical case - first argument is valid
    def test_df_to_pandas(self):
        data = StringIO("# Accession:\n"
                        "# Downloaded (UTC):\n"
                        "# Licence:\n"
                        "# Licence URL:\n"
                        "accession,hgvs_nt,hgvs_pro,score\n"
                        "urn:mavedb:00000011-a-1#27,c.[4=;5=;6=],p.Phe2=,0.0")
        results = df_to_pandas(data, True)
        self.assertTrue(results[0].iat[0, 2] == "p.Phe2=")

    # scenario 2 - call function on invalid filename

    # edge case - invalid first argument
    def test_df_to_pandas_invalid_first_arg(self):
        with self.assertRaises(ValueError):
            df_to_pandas("notarealfile")

    # scenario 3 - call function without arguments

    # edge case - missing first argument
    def test_df_to_pandas_missing_arguments(self):
        with self.assertRaises(TypeError):
            df_to_pandas()
