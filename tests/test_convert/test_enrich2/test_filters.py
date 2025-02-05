import unittest

import numpy as np
import pandas as pd

from mavetools.convert.enrich2 import constants, filters


class TestDropNaColumns(unittest.TestCase):
    def test_drops_null_nt_column(self):
        df = pd.DataFrame(
            {
                constants.nt_variant_col: [None, None],
                constants.pro_variant_col: ["p.G4L", "p.G4L"],
            }
        )
        filters.drop_na_columns(df)
        self.assertNotIn(constants.nt_variant_col, df)

    def test_drops_null_pro_column(self):
        df = pd.DataFrame(
            {
                constants.pro_variant_col: [None, None],
                constants.nt_variant_col: ["c.100A>G", "c.101A>G"],
            }
        )
        filters.drop_na_columns(df)
        self.assertNotIn(constants.pro_variant_col, df)

    def test_does_not_drop_partial_pro_column(self):
        df = pd.DataFrame(
            {
                constants.pro_variant_col: [None, "pG4L"],
                constants.nt_variant_col: ["c.100A>G", "c.101A>G"],
            }
        )
        filters.drop_na_columns(df)
        self.assertIn(constants.pro_variant_col, df)

    def test_does_not_drop_nt_column_when_defined(self):
        df = pd.DataFrame(
            {
                constants.pro_variant_col: [None, "pG4L"],
                constants.nt_variant_col: ["c.100A>G", "c.101A>G"],
            }
        )
        filters.drop_na_columns(df)
        self.assertIn(constants.nt_variant_col, df)

    def test_drops_null_column(self):
        df = pd.DataFrame({"A": [None, np.NaN]})
        filters.drop_na_columns(df)
        self.assertNotIn("A", df)

    def test_does_not_drop_column_containing_non_null_values(self):
        df = pd.DataFrame({"A": [None, "pG4L"]})
        filters.drop_na_columns(df)
        self.assertIn("A", df)


class TestDropNaRows(unittest.TestCase):
    def test_drops_null_row(self):
        df = pd.DataFrame({"A": [None], "B": [np.NaN]})
        filters.drop_na_rows(df)
        self.assertEqual(len(df), 0)

    def test_does_not_drop_row_containing_non_null_values(self):
        df = pd.DataFrame({"A": [None], "B": [0.0]})
        filters.drop_na_rows(df)
        self.assertEqual(len(df), 1)


if __name__ == "__main__":
    unittest.main()
