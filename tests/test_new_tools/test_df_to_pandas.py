from unittest import TestCase
import pandas as pd
from mavetools.new_tools.df_to_pandas import df_to_pandas


class Test(TestCase):

    # scenario 1 - call function on urn mavedb 00000001-a-1_scores

    # typical case - first argument is valid
    def test_df_to_pandas(self):
        results = df_to_pandas("urn mavedb 00000001-a-1_scores.csv")
        self.assertTrue(results[0].iat[0, 2] == "p.Pro73Gln")

    # typical case - first argument is valid, second argument passed
    def test_df_to_pandas_drop_accession(self):
        results = df_to_pandas("urn mavedb 00000001-a-1_scores.csv", True)
        self.assertTrue(results[0].iat[0, 1] == "p.Pro73Gln")

    # edge case - invalid second argument
    def test_df_to_pandas_invalid_second_arg(self):
        with self.assertRaises(TypeError):
            df_to_pandas("urn mavedb 00000001-a-1_scores.csv", 4)

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
