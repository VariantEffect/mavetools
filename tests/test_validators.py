import unittest

import pandas as pd

from mavedbconvert import validators, constants, exceptions


class TestHGVSPatternsBackend(unittest.TestCase):
    def setUp(self):
        self.backend = validators.HGVSPatternsBackend()

    def test_validate_hgvs_raise_HGVSValidationError(self):
        with self.assertRaises(exceptions.HGVSValidationError):
            self.backend.validate("p.1102A>G")
        with self.assertRaises(exceptions.HGVSValidationError):
            self.backend.validate("x.102A>G")

    def test_validate_passes_on_special(self):
        self.backend.validate(constants.enrich2_wildtype)
        self.backend.validate(constants.enrich2_synonymous)

    def test_returns_str_variant(self):
        self.assertIsInstance(self.backend.validate("c.1A>G"), str)


class TestValidateHGVS(unittest.TestCase):
    def test_uses_patterns_backend_as_default(self):
        result = validators.validate_variants(["c.[1A>G;2A>G]"], n_jobs=2, verbose=0)
        self.assertIsInstance(result[0], str)

    def test_can_specify_backend(self):
        backend = validators.HGVSPatternsBackend()
        result = validators.validate_variants(
            ["c.[1A>G;2A>G]"], n_jobs=2, verbose=0, validation_backend=backend
        )
        self.assertIsInstance(result[0], str)


class TestDfValidators(unittest.TestCase):
    def test_validate_column_raise_keyerror_column_not_exist(self):
        df = pd.DataFrame({"a": [1]})
        with self.assertRaises(KeyError):
            validators.validate_has_column(df, "b")

    def test_validate_column_passes_when_column_exists(self):
        df = pd.DataFrame({"a": [1]})
        validators.validate_has_column(df, "a")

    def test_error_some_values_non_numeric(self):
        df = pd.DataFrame({"A": ["a", 1, 2]})
        with self.assertRaises(TypeError):
            validators.validate_columns_are_numeric(df)

    def test_pass_all_numeric(self):
        df = pd.DataFrame({"A": [1, 2, 1.0]})
        validators.validate_columns_are_numeric(df)


class TestHGVSValidators(unittest.TestCase):
    def test_validate_hgvs_uniqueness(self):
        df = pd.DataFrame({constants.nt_variant_col: ["a", "b"]})
        validators.validate_hgvs_uniqueness(df, constants.nt_variant_col)  # Should pass

        df = pd.DataFrame({constants.nt_variant_col: ["a", "b", "a"]})
        with self.assertRaises(ValueError):
            validators.validate_hgvs_uniqueness(df, constants.nt_variant_col)

        # test multi-variant formatting
        df = pd.DataFrame({constants.nt_variant_col: list("abcdefg" * 2)})
        with self.assertRaises(ValueError) as cm:
            validators.validate_hgvs_uniqueness(df, constants.nt_variant_col)
        self.assertTrue(str(cm.exception).endswith(", ..."))

    def test_validate_hgvs_uniqueness_bad_column(self):
        df = pd.DataFrame({constants.nt_variant_col: ["a", "b", "a"]})
        with self.assertRaises(KeyError):
            validators.validate_hgvs_uniqueness(df, constants.pro_variant_col)

    def test_validate_hgvs_uniqueness_ignores_none(self):
        df = pd.DataFrame({constants.nt_variant_col: ["a", "b", None, None]})
        validators.validate_hgvs_uniqueness(df, constants.nt_variant_col)  # Should pass


class TestMaveDBCompliance(unittest.TestCase):
    def test_error_primary_column_contains_null(self):
        df = pd.DataFrame(
            {
                constants.nt_variant_col: ["c.100A>G", None],
                constants.pro_variant_col: ["p.G4L", "p.G5L"],
            }
        )
        with self.assertRaises(ValueError):
            validators.validate_mavedb_compliance(df, df_type=None)

    def test_error_primary_column_as_pro_contains_null(self):
        df = pd.DataFrame(
            {
                constants.nt_variant_col: [None, None],
                constants.pro_variant_col: ["p.G4L", None],
            }
        )
        with self.assertRaises(ValueError):
            validators.validate_mavedb_compliance(df, df_type=None)

    def test_pass_coding_(self):
        df = pd.DataFrame(
            {
                constants.nt_variant_col: ["c.100A>G", "c.101A>G"],
                constants.pro_variant_col: ["p.G4L", "p.G5L"],
            }
        )
        validators.validate_mavedb_compliance(df, df_type=None)

        df = pd.DataFrame(
            {
                constants.nt_variant_col: ["n.100A>G", "n.101A>G"],
                constants.pro_variant_col: [None, None],
            }
        )
        validators.validate_mavedb_compliance(df, df_type=None)

    def test_error_missing_nt_pro_columns(self):
        df = pd.DataFrame({"A": ["c.100A>G", "c.101A>G"], "B": [None, None]})
        with self.assertRaises(ValueError):
            validators.validate_mavedb_compliance(df, df_type=None)

    def test_error_neither_column_defines_variants(self):
        df = pd.DataFrame(
            {
                constants.nt_variant_col: [None, None],
                constants.pro_variant_col: [None, None],
            }
        )
        with self.assertRaises(ValueError):
            validators.validate_mavedb_compliance(df, df_type=None)

    def test_allows_duplicates_in_pro_col(self):
        df = pd.DataFrame(
            {
                constants.nt_variant_col: [None, None],
                constants.pro_variant_col: ["p.G4L", "p.G4L"],
            }
        )
        validators.validate_mavedb_compliance(df, df_type=None)  # passes

    def test_error_duplicates_in_nt_col(self):
        df = pd.DataFrame(
            {
                constants.nt_variant_col: ["c.100A>G", "c.100A>G"],
                constants.pro_variant_col: ["p.G4L", "p.G4L"],
            }
        )
        with self.assertRaises(ValueError):
            validators.validate_mavedb_compliance(df, df_type=None)

    def test_keyerror_missing_score_column_df_type_is_scores(self):
        df = pd.DataFrame(
            {
                constants.pro_variant_col: [None, "pG4L"],
                constants.nt_variant_col: ["c.100A>G", "c.101A>G"],
            }
        )
        with self.assertRaises(KeyError):
            validators.validate_mavedb_compliance(df, df_type=constants.score_type)


class TestValidateSameVariants(unittest.TestCase):
    def test_ve_counts_defines_different_nt_variants(self):
        scores = pd.DataFrame(
            {
                constants.nt_variant_col: ["c.1A>G"],
                constants.pro_variant_col: ["p.Leu5Glu"],
            }
        )
        counts = pd.DataFrame(
            {
                constants.nt_variant_col: ["c.2A>G"],
                constants.pro_variant_col: ["p.Leu5Glu"],
            }
        )
        with self.assertRaises(AssertionError):
            validators.validate_datasets_define_same_variants(scores, counts)

        scores = pd.DataFrame({constants.nt_variant_col: ["n.1A>G"]})
        counts = pd.DataFrame({constants.nt_variant_col: ["n.2A>G"]})
        with self.assertRaises(AssertionError):
            validators.validate_datasets_define_same_variants(scores, counts)

    def test_ve_counts_defines_different_pro_variants(self):
        scores = pd.DataFrame(
            {
                constants.nt_variant_col: ["c.1A>G"],
                constants.pro_variant_col: ["p.Leu5Glu"],
            }
        )
        counts = pd.DataFrame(
            {
                constants.nt_variant_col: ["c.1A>G"],
                constants.pro_variant_col: ["p.Leu75Glu"],
            }
        )
        with self.assertRaises(AssertionError):
            validators.validate_datasets_define_same_variants(scores, counts)

        scores = pd.DataFrame({constants.pro_variant_col: ["p.Leu5Glu"]})
        counts = pd.DataFrame({constants.pro_variant_col: ["p.Leu75Glu"]})
        with self.assertRaises(AssertionError):
            validators.validate_datasets_define_same_variants(scores, counts)

    def test_passes_when_same_variants_defined(self):
        scores = pd.DataFrame(
            {
                constants.nt_variant_col: ["c.1A>G"],
                constants.pro_variant_col: ["p.Leu5Glu"],
            }
        )
        counts = pd.DataFrame(
            {
                constants.nt_variant_col: ["c.1A>G"],
                constants.pro_variant_col: ["p.Leu5Glu"],
            }
        )
        validators.validate_datasets_define_same_variants(scores, counts)

        scores = pd.DataFrame({constants.nt_variant_col: ["n.1A>G"]})
        counts = pd.DataFrame({constants.nt_variant_col: ["n.1A>G"]})
        validators.validate_datasets_define_same_variants(scores, counts)

        scores = pd.DataFrame({constants.pro_variant_col: ["p.Leu5Glu"]})
        counts = pd.DataFrame({constants.pro_variant_col: ["p.Leu5Glu"]})
        validators.validate_datasets_define_same_variants(scores, counts)

    def test_error_dfs_define_different_hgvs_columns(self):
        scores = pd.DataFrame({constants.nt_variant_col: ["c.1A>G"]})
        counts = pd.DataFrame({constants.pro_variant_col: ["p.Leu75Glu"]})
        with self.assertRaises(AssertionError):
            validators.validate_datasets_define_same_variants(scores, counts)


if __name__ == "__main__":
    unittest.main()