import os
import unittest
from unittest.mock import patch

import re

from mavehgvs.patterns import dna, protein

import numpy as np
import pandas as pd
from pandas.testing import assert_index_equal
from fqfa.constants.translation.table import CODON_TABLE
from fqfa.constants.iupac.protein import AA_CODES

from mavetools.mavedbconvert import enrich2, constants, exceptions

from tests import ProgramTestCase


# Utility tests
# --------------------------------------------------------------------------- #
class TestGetCountDataFrames(ProgramTestCase):
    """
    Test method get_count_dataframes checking if conditions are correctly
    parsed.
    """

    def setUp(self):
        super().setUp()
        self.path = os.path.join(self.data_dir, "enrich2", "test_store.h5")
        self.store = pd.HDFStore(self.path, "w")
        index = pd.MultiIndex.from_product(
            [["c1", "c2"], ["rep1", "rep2"], ["t0", "t1"]],
            names=["condition", "selection", "timepoint"],
        )
        scores_hgvs = ["c.1A>G", "c.3A>G"]
        counts_hgvs = ["c.1A>G", "c.2A>G"]
        self.store["/main/variants/scores/"] = pd.DataFrame(
            np.random.randn(len(scores_hgvs), len(index)),
            index=scores_hgvs,
            columns=index,
        )
        self.store["/main/variants/counts/"] = pd.DataFrame(
            np.random.randn(len(scores_hgvs), len(index)),
            index=counts_hgvs,
            columns=index,
        )

    def tearDown(self):
        self.store.close()
        if os.path.isfile(self.path):
            os.unlink(self.path)

    def test_column_names_combine_selection_and_timepoint(self):
        cnd_df = enrich2.get_count_dataframe_by_condition(self.store, cnd="c1")
        self.assertListEqual(
            list(cnd_df.columns), ["rep1_t0", "rep1_t1", "rep2_t0", "rep2_t1"]
        )

    def test_index_of_dfs_match_index_of_scores(self):
        cnd_df = enrich2.get_count_dataframe_by_condition(self.store, cnd="c1")
        assert_index_equal(self.store["/main/variants/scores/"].index, cnd_df.index)

    def test_row_filled_with_nans_filtered_index_not_in_counts(self):
        cnd_df = enrich2.get_count_dataframe_by_condition(self.store, cnd="c1")
        self.assertTrue(np.all(cnd_df.loc["c.3A>G", :].isnull()))

    def test_returns_empty_when_missing_scores_key(self):
        self.store.remove("/main/variants/scores/")
        cnd_df = enrich2.get_count_dataframe_by_condition(self.store, cnd="c1")
        self.assertIsNone(cnd_df)

    def test_returns_empty_when_missing_counts_key(self):
        self.store.remove("/main/variants/counts/")
        cnd_df = enrich2.get_count_dataframe_by_condition(self.store, cnd="c1")
        self.assertIsNone(cnd_df)


class TestFlattenColumnNames(unittest.TestCase):
    def setUp(self):
        index = pd.MultiIndex.from_product(
            [["c1", "c2"], ["rep1", "rep2"], ["t0", "t1"]],
            names=["condition", "selection", "timepoint"],
        )
        scores_hgvs = ["c.1A>G", "c.3A>G"]
        self.df = pd.DataFrame(
            np.random.randn(len(scores_hgvs), len(index)),
            index=scores_hgvs,
            columns=index,
        )

    def test_column_names_combine_columns_using_ordering(self):
        cnames = enrich2.flatten_column_names(
            self.df.loc[:, pd.IndexSlice["c1", :, :]].columns, ordering=(2, 1)
        )
        self.assertListEqual(cnames, ["t0_rep1", "t1_rep1", "t0_rep2", "t1_rep2"])


class TestReplicateScoreDataFrames(ProgramTestCase):
    """
    Test method get_replicate_score_dataframes checking if conditions are
    correctly parsed.
    """

    def setUp(self):
        super().setUp()
        self.path = os.path.join(self.data_dir, "enrich2", "test_store.h5")
        self.store = pd.HDFStore(self.path, "w")

        shared_index = pd.MultiIndex.from_product(
            [["c1", "c2"], ["rep1", "rep2"], ["SE", "score"]],
            names=["condition", "selection", "value"],
        )
        index = pd.MultiIndex.from_product(
            [["c1", "c2"], ["SE", "epsilon", "score"]], names=["condition", "value"]
        )

        hgvs = ["c.1A>G", "c.2A>G"]
        self.store["/main/variants/scores/"] = pd.DataFrame(
            np.random.randn(len(hgvs), len(index)), index=hgvs, columns=index
        )
        self.store["/main/variants/scores_shared/"] = pd.DataFrame(
            np.random.randn(len(hgvs), len(shared_index)),
            index=hgvs,
            columns=shared_index,
        )

    def tearDown(self):
        super().tearDown()
        self.store.close()
        if os.path.isfile(self.path):
            os.unlink(self.path)

    def test_conditions_are_dictionary_keys(self):
        cnd_dfs = enrich2.get_replicate_score_dataframes(self.store)
        self.assertIn("c1", cnd_dfs)
        self.assertIn("c2", cnd_dfs)

    def test_returns_empty_when_missing_scores_key(self):
        self.store.remove("/main/variants/scores")
        cnd_dfs = enrich2.get_replicate_score_dataframes(self.store)
        self.assertDictEqual(cnd_dfs, {})

    def test_returns_empty_when_missing_shared_scores_key(self):
        self.store.remove("/main/variants/scores_shared")
        cnd_dfs = enrich2.get_replicate_score_dataframes(self.store)
        self.assertDictEqual(cnd_dfs, {})

    def test_adds_rep_id_to_score_and_SE(self):
        cnd_dfs = enrich2.get_replicate_score_dataframes(self.store)
        for cnd, df in cnd_dfs.items():
            for c_name in df.columns:
                if c_name not in ("score", "SE", "epsilon"):
                    self.assertIn("rep", c_name.lower())

    def test_assertion_error_scores_shared_scores_different_index(self):
        hgvs = ["c.1A>G", "c.3A>G"]
        shared_index = pd.MultiIndex.from_product(
            [["c1", "c2"], ["rep1", "rep2"], ["SE", "score"]],
            names=["condition", "selection", "value"],
        )
        self.store["/main/variants/scores_shared/"] = pd.DataFrame(
            np.random.randn(len(hgvs), len(shared_index)),
            index=hgvs,
            columns=shared_index,
        )
        with self.assertRaises(AssertionError):
            enrich2.get_replicate_score_dataframes(self.store)


class TestDropNull(unittest.TestCase):
    def test_calls_drop_na_rows_from_scores_inplace(self):
        df = pd.DataFrame({"A": [None, 1]})
        enrich2.drop_null(df)
        assert_index_equal(df.index, pd.Index([1]))

    def test_calls_drop_na_cols_from_scores_inplace(self):
        df = pd.DataFrame({"A": [1, 2], "B": [None, None]})
        enrich2.drop_null(df)
        self.assertNotIn("B", df)

    def test_assertion_error_if_counts_scores_indices_do_not_match(self):
        df1 = pd.DataFrame({"A": [1, 2]}, index=["a", "b"])
        df2 = pd.DataFrame({"A": [1, 2]}, index=["c", "b"])
        with self.assertRaises(AssertionError):
            enrich2.drop_null(df1, df2)

    def test_assertion_error_if_dfs_define_different_variants(self):
        df1 = pd.DataFrame(
            {
                constants.nt_variant_col: ["c.1A>G"],
                constants.pro_variant_col: ["p.G4L"],
            },
            index=["c.1A>G"],
        )
        df2 = pd.DataFrame(
            {
                constants.nt_variant_col: ["c.2A>G"],
                constants.pro_variant_col: ["p.G4L"],
            },
            index=["c.1A>G"],
        )
        with self.assertRaises(AssertionError):
            enrich2.drop_null(df1, df2)

    def test_na_rows_dropped_from_scores_counts_after_join(self):
        df1 = pd.DataFrame(
            {
                constants.nt_variant_col: ["c.1A>G", "c.2A>G"],
                constants.pro_variant_col: ["p.G4L", "p.G5L"],
                "score": [1, None],
            },
            index=["c.1A>G", "c.2A>G"],
        )
        df2 = pd.DataFrame(
            {
                constants.nt_variant_col: ["c.1A>G", "c.2A>G"],
                constants.pro_variant_col: ["p.G4L", "p.G5L"],
                "count": [10, None],
            },
            index=["c.1A>G", "c.2A>G"],
        )
        scores, counts = enrich2.drop_null(df1, df2)
        self.assertNotIn("c.2A>G", scores.index.values)
        self.assertNotIn("c.2A>G", counts.index.values)

    def test_na_cols_dropped_from_scores_counts_after_join(self):
        df1 = pd.DataFrame(
            {
                constants.nt_variant_col: ["c.1A>G"],
                constants.pro_variant_col: [None],
                "score": [1],
            },
            index=["c.1A>G"],
        )
        df2 = pd.DataFrame(
            {
                constants.nt_variant_col: ["c.1A>G"],
                constants.pro_variant_col: [None],
                "count": [10],
            },
            index=["c.1A>G"],
        )
        scores, counts = enrich2.drop_null(df1, df2)
        self.assertNotIn(constants.pro_variant_col, scores.columns)
        self.assertNotIn(constants.pro_variant_col, counts.columns)

    def test_scores_and_counts_columns_separated_after_join(self):
        df1 = pd.DataFrame(
            {
                constants.nt_variant_col: ["c.1A>G", "c.2A>G"],
                constants.pro_variant_col: ["p.G4L", "p.G5L"],
                "score": [1, None],
            },
            index=["c.1A>G", "c.2A>G"],
        )
        df2 = pd.DataFrame(
            {
                constants.nt_variant_col: ["c.1A>G", "c.2A>G"],
                constants.pro_variant_col: ["p.G4L", "p.G5L"],
                "count": [10, None],
            },
            index=["c.1A>G", "c.2A>G"],
        )
        scores, counts = enrich2.drop_null(df1, df2)
        self.assertListEqual(
            list(scores.columns),
            [constants.nt_variant_col, constants.pro_variant_col, "score"],
        )
        self.assertListEqual(
            list(counts.columns),
            [constants.nt_variant_col, constants.pro_variant_col, "count"],
        )


# HD5/Row parsing tests
# --------------------------------------------------------------------------- #
class TestEnrich2ConvertH5Filepath(ProgramTestCase):
    def setUp(self):
        super().setUp()
        self.path = os.path.join(self.data_dir, "enrich2", "enrich2.h5")
        self.enrich2 = enrich2.Enrich2(self.path, wt_sequence="AAA")

    def test_replaces_underscore_with_spaces(self):
        res = self.enrich2.convert_h5_filepath("base", "syn vars", "scores", "c1")
        self.assertIn("syn_vars", res)

    def test_concats_basename_elem_type_then_cnd_and_csv_ext(self):
        res = self.enrich2.convert_h5_filepath("base", "syn vars", "scores", "c1")
        expected = "mavedb_base_syn_vars_scores_c1.csv"
        self.assertEqual(res, os.path.join(self.enrich2.output_directory, expected))


class TestEnrich2ConvertH5Df(ProgramTestCase):
    def setUp(self):
        super().setUp()
        self.path = os.path.join(self.data_dir, "enrich2", "enrich2.h5")
        self.enrich2 = enrich2.Enrich2(self.path, wt_sequence="AAA")

    def test_doesnt_open_invalid_rows_file_if_there_are_no_invalid_rows(self):
        self.path = os.path.join(self.data_dir, "enrich2", "enrich2.tsv")
        self.enrich2 = enrich2.Enrich2(self.path, wt_sequence="AAA")
        invalid_rows_path = os.path.join(
            os.path.dirname(self.path), "enrich2_invalid_rows.csv"
        )

        df = pd.DataFrame(data={"score": [1]}, index=["c.1A>G (p.Lys1Val)"])
        self.enrich2.convert_h5_df(
            df=df, element=constants.variants_table, df_type=constants.score_type
        )
        self.assertFalse(os.path.isfile(invalid_rows_path))

    def test_drops_non_numeric_columns(self):
        df = pd.DataFrame(data={"score": [1], "B": ["a"]}, index=["c.1A>G (p.Lys1Val)"])
        result = self.enrich2.convert_h5_df(
            df=df, element=constants.variants_table, df_type=constants.score_type
        )
        self.assertNotIn("B", result)

    def test_type_casts_numeric_to_int_and_float(self):
        df = pd.DataFrame(data={"score": [1], "B": [1.2]}, index=["c.1A>G (p.Lys1Val)"])
        result = self.enrich2.convert_h5_df(
            df=df, element=constants.variants_table, df_type=constants.score_type
        )
        self.assertTrue(np.issubdtype(result["score"].values[0], np.signedinteger))
        self.assertTrue(np.issubdtype(result["B"].values[0], np.floating))

    def test_sets_index_as_input_index(self):
        df = pd.DataFrame({"score": [1], "B": ["a"]}, index=["c.1A>T (p.Lys1Val)"])
        result = self.enrich2.convert_h5_df(
            df=df, element=constants.variants_table, df_type=constants.score_type
        )
        assert_index_equal(result.index, df.index)

    def test_opens_invalid_rows_file_for_invalid_rows(self):
        self.path = os.path.join(self.data_dir, "enrich2", "enrich2.tsv")
        self.enrich2 = enrich2.Enrich2(self.path, wt_sequence="AAA")
        df = pd.DataFrame(data={"score": [1], "B": ["a"]}, index=["c.1T>G (p.Lys1Val)"])
        with self.assertRaises(ValueError):
            self.enrich2.convert_h5_df(
                df=df, element=constants.variants_table, df_type=constants.score_type
            )

        invalid_rows_path = os.path.join(
            os.path.dirname(self.path), "enrich2_invalid_rows.csv"
        )

        self.assertTrue(os.path.isfile(invalid_rows_path))

    def test_opens_invalid_rows_file_for_invalid_rows_by_condition(self):
        self.path = os.path.join(self.data_dir, "enrich2", "enrich2.tsv")
        self.enrich2 = enrich2.Enrich2(self.path, wt_sequence="AAA")
        df = pd.DataFrame(data={"score": [1], "B": ["a"]}, index=["c.1T>G (p.Lys1Val)"])
        with self.assertRaises(ValueError):
            self.enrich2.convert_h5_df(
                df=df,
                element=constants.variants_table,
                df_type=constants.score_type,
                cnd="c1",
            )

        invalid_rows_path = os.path.join(
            os.path.dirname(self.path),
            "mavedb_enrich2_variants_counts_c1_invalid_rows.csv",
        )

        self.assertTrue(os.path.isfile(invalid_rows_path))

    def test_invalid_rows_file_contains_error_description(self):
        self.path = os.path.join(self.data_dir, "enrich2", "enrich2.tsv")
        self.enrich2 = enrich2.Enrich2(self.path, wt_sequence="AAA")
        invalid_rows_path = os.path.join(
            os.path.dirname(self.path), "enrich2_invalid_rows.csv"
        )

        df = pd.DataFrame(
            data={"score": [1.1, 1.2]},
            index=["c.1A>T (p.Lys1Val)", "c.1T>G (p.Lys1Val)"],
        )

        self.enrich2.convert_h5_df(df=df, df_type=constants.score_type, element=None)
        self.assertTrue(os.path.isfile(invalid_rows_path))

        invalid = pd.read_csv(invalid_rows_path, sep=",", index_col=0)
        self.assertEqual(len(invalid), 1)
        self.assertEqual(invalid.index[0], "c.1T>G (p.Lys1Val)")
        self.assertIn("error_description", invalid.columns)


class TestEnrich2LoadInput(ProgramTestCase):
    def test_error_file_not_h5_or_tsv(self):
        path = os.path.join(self.data_dir, "empiric", "empiric.xlsx")
        p = enrich2.Enrich2(path, wt_sequence="AAA")
        with self.assertRaises(TypeError):
            p.load_input_file()

    def test_scores_tsv_missing_score_column(self):
        path = os.path.join(self.data_dir, "enrich2", "enrich2.tsv")
        p = enrich2.Enrich2(
            path,
            wt_sequence="AAA",
            score_column="scores",
            hgvs_column="sequence",
            input_type=constants.score_type,
        )
        with self.assertRaises(KeyError):
            p.load_input_file()

    def test_input_type_counts_doesnt_raise_keyerror(self):
        path = os.path.join(self.data_dir, "enrich2", "enrich2.tsv")
        p = enrich2.Enrich2(
            path,
            wt_sequence="AAA",
            hgvs_column="sequence",
            input_type=constants.count_type,
        )
        p.load_input_file()

    def test_scores_tsv_missing_hgvs_column(self):
        path = os.path.join(self.data_dir, "enrich2", "enrich2.tsv")
        p = enrich2.Enrich2(path, wt_sequence="AAA", hgvs_column="hgvs")
        with self.assertRaises(KeyError):
            p.load_input_file()


class TestEnrich2ParseRow(ProgramTestCase):
    def setUp(self):
        super().setUp()
        self.path = os.path.join(self.data_dir, "enrich2", "dummy.h5")
        self.enrich2 = enrich2.Enrich2(self.path, wt_sequence="ACT")

    def test_invalid_variant(self):
        with self.assertRaises(exceptions.InvalidVariantType):
            self.enrich2.parse_row(("c.1_2del", None))
        with self.assertRaises(exceptions.InvalidVariantType):
            self.enrich2.parse_row(("b.1A>G", None))

    def test_special_nt_variants_are_singletons(self):
        with self.assertRaises(ValueError):
            self.enrich2.parse_nucleotide_variant("_wt, n.2C>G")
        with self.assertRaises(ValueError):
            self.enrich2.parse_nucleotide_variant("n.2C>G, _wt")

    def test_special_pro_variants_are_singletons(self):
        with self.assertRaises(ValueError):
            self.enrich2.parse_protein_variant("_sy, p.Thr1Gly")
        with self.assertRaises(ValueError):
            self.enrich2.parse_protein_variant("p.Thr1Gly, _sy")

    def test_infers_dna(self):
        # test tuples
        for prefix in "cngm":
            variant = "{0}.1A>G, {0}.2C>G".format(prefix)
            expected = "{0}.[1A>G;2C>G]".format(prefix), None
            self.assertEqual(expected, self.enrich2.parse_row((variant, None)))

        # test strings
        for prefix in "cngm":
            variant = "{0}.1A>G, {0}.2C>G".format(prefix)
            expected = "{0}.[1A>G;2C>G]".format(prefix), None
            self.assertEqual(expected, self.enrich2.parse_row(variant))

    def test_invalid_multiprefix(self):
        with self.assertRaises(ValueError):
            self.enrich2.parse_row("c.1A>G, n.2C>G")

    def test_infers_protein(self):
        # test tuples
        variant = "p.Thr1=, p.Thr1Gly"
        expected = (None, "p.[Thr1=;Thr1Gly]")
        self.assertEqual(expected, self.enrich2.parse_row((variant, None)))

        # test strings
        variant = "p.Thr1=, p.Thr1Gly"
        expected = (None, "p.[Thr1=;Thr1Gly]")
        self.assertEqual(expected, self.enrich2.parse_row(variant))

    def test_nt_variant_is_none_special_variant_is_from_synonymous_table(self):
        self.assertEqual(
            (None, constants.enrich2_synonymous),
            self.enrich2.parse_row(
                (constants.enrich2_synonymous, constants.synonymous_table)
            ),
        )

    @patch("mavetools.mavedbconvert.enrich2.apply_offset", return_value="c.3T>C (p.Thr1=)")
    def test_calls_apply_offset_to_variant(self, patch):
        variant = "c.3T>C (p.=)"
        self.enrich2.parse_row((variant, None))
        patch.assert_called()

    def test_delegate_to_multi(self):
        variant = "c.3T>C (p.Thr1=)"
        expected = ("c.3T>C", "p.Thr1=")
        self.assertEqual(expected, self.enrich2.parse_row((variant, None)))

    def test_returns_special_variant_as_tuple_non_synonymous_table(self):
        self.assertEqual(self.enrich2.parse_row(("_wt", None)), ("_wt", "_wt"))
        self.assertEqual(self.enrich2.parse_row(("_sy", None)), ("_sy", "_sy"))

    def test_strips_whitespace(self):
        self.assertEqual(self.enrich2.parse_row((" c.1A>G ", None)), ("c.1A>G", None))

    # TODO: Uncomment if normalizing variants.
    # def test_converts_X_to_N(self):
    #     self.assertEqual(
    #         self.enrich2.parse_row(("c.1A>X", None)), ("c.1A>N", None)
    #     )

    # TODO: Uncomment if normalizing variants.
    # def test_converts_qmarks_to_Xaa_to_single_q(self):
    #     self.assertEqual(
    #         self.enrich2.parse_row(("p.Thr1???", None)), (None, "p.Thr1Xaa")
    #     )
    #     self.assertEqual(
    #         self.enrich2.parse_row(("p.T1?", None)), (None, "p.Thr1Xaa")
    #     )


# Protein parsing tests
# --------------------------------------------------------------------------- #
class TestProteinHGVSParsing(ProgramTestCase):
    def setUp(self):
        super().setUp()
        self.path = os.path.join(self.data_dir, "enrich2", "dummy.h5")
        self.enrich2 = enrich2.Enrich2(self.path, wt_sequence="AAA")

    def test_combines_string_into_multi_variant_syntax(self):
        result = self.enrich2.parse_protein_variant("p.Leu5Gly,p.Leu6Gly")
        self.assertEqual(result, "p.[Leu5Gly;Leu6Gly]")

    def test_combines_list_into_multi_variant_syntax(self):
        result = self.enrich2.parse_protein_variant(["p.Leu5Gly", "p.Leu6Gly"])
        self.assertEqual(result, "p.[Leu5Gly;Leu6Gly]")

    def test_does_not_squash_single_variant(self):
        result = self.enrich2.parse_protein_variant("p.Leu5Gly")
        self.assertEqual(result, "p.Leu5Gly")

    def test_passes_on_sy_or_wt(self):
        self.assertEqual(self.enrich2.parse_protein_variant("_wt"), "_wt")
        self.assertEqual(self.enrich2.parse_protein_variant("_sy"), "_sy")

    def test_valueerr_not_a_valid_protein_syntax(self):
        with self.assertRaises(ValueError):
            self.enrich2.parse_protein_variant("c.101A>G")
        with self.assertRaises(ValueError):
            self.enrich2.parse_protein_variant("p.101A>G")
        with self.assertRaises(ValueError):
            self.enrich2.parse_protein_variant("random")
        with self.assertRaises(ValueError):
            self.enrich2.parse_protein_variant("(random)")

    def test_removes_brackets(self):
        result = self.enrich2.parse_protein_variant("(p.Leu4Gly),(p.Leu5Gly)")
        self.assertEqual(result, "p.[Leu4Gly;Leu5Gly]")

        result = self.enrich2.parse_protein_variant("(p.Leu4Gly)")
        self.assertEqual(result, "p.Leu4Gly")

        result = self.enrich2.parse_protein_variant(["(p.Leu4Gly)"])
        self.assertEqual(result, "p.Leu4Gly")

    def test_strips_ws(self):
        result = self.enrich2.parse_protein_variant(" p.Gly5Leu ")
        self.assertEqual(result, "p.Gly5Leu")
        result = self.enrich2.parse_protein_variant(" p.Leu4Gly, p.Leu5Gly ")
        self.assertEqual(result, "p.[Leu4Gly;Leu5Gly]")

    def test_removes_duplicates(self):
        result = self.enrich2.parse_protein_variant("p.Leu5Gly,p.Leu5Gly")
        self.assertEqual(result, "p.Leu5Gly")

    def test_maintains_ordering(self):
        result = self.enrich2.parse_protein_variant("p.Leu5Gly,p.Leu4Gly")
        self.assertEqual(result, "p.[Leu5Gly;Leu4Gly]")


# Nucleotide parsing tests
# --------------------------------------------------------------------------- #
class TestNucleotideHGVSParing(ProgramTestCase):
    def setUp(self):
        super().setUp()
        self.path = os.path.join(self.data_dir, "enrich2", "dummy.h5")
        self.enrich2 = enrich2.Enrich2(self.path, wt_sequence="AAA")

    def test_parses_non_coding_nt_variants_into_multi_variant(self):
        nt = self.enrich2.parse_nucleotide_variant(
            "c.-455T>A, c.-122A>T, c.-101A>T, c.-42T>A"
        )
        self.assertEqual(nt, "c.[-455T>A;-122A>T;-101A>T;-42T>A]")

    def test_combines_into_multi_variant_syntax(self):
        result = self.enrich2.parse_nucleotide_variant("c.2A>G,c.1A>G")
        self.assertEqual(result, "c.[2A>G;1A>G]")

    def test_combines_list_into_multi_variant_syntax(self):
        result = self.enrich2.parse_nucleotide_variant(["c.2A>G", "c.1A>G"])
        self.assertEqual(result, "c.[2A>G;1A>G]")

    def test_does_not_squash_single_variant(self):
        result = self.enrich2.parse_nucleotide_variant("c.2A>G")
        self.assertEqual(result, "c.2A>G")

    def test_passes_on_sy_or_wt(self):
        self.assertEqual(self.enrich2.parse_nucleotide_variant("_wt"), "_wt")
        self.assertEqual(self.enrich2.parse_nucleotide_variant("_sy"), "_sy")

    def test_valueerr_not_a_valid_protein_syntax(self):
        with self.assertRaises(ValueError):
            self.enrich2.parse_nucleotide_variant("p.101A>G")
        with self.assertRaises(ValueError):
            self.enrich2.parse_nucleotide_variant("p.Leu5Gly")
        with self.assertRaises(ValueError):
            self.enrich2.parse_nucleotide_variant("random")
        with self.assertRaises(ValueError):
            self.enrich2.parse_nucleotide_variant("()")
        with self.assertRaises(ValueError):
            self.enrich2.parse_nucleotide_variant("c.[2A>G;1A>G]")

    def test_valueerr_multi_prefix_types(self):
        with self.assertRaises(ValueError):
            self.enrich2.parse_nucleotide_variant("c.1A>G;n.2A>G")

    def test_strips_ws(self):
        result = self.enrich2.parse_nucleotide_variant(" c.101A>G ")
        self.assertEqual(result, "c.101A>G")
        result = self.enrich2.parse_nucleotide_variant(" c.2A>G, c.2A>G ")
        self.assertEqual(result, "c.[2A>G;2A>G]")


# Mixed parsing tests
# --------------------------------------------------------------------------- #
class TestEnrich2MixedHGVSParsing(ProgramTestCase):
    def setUp(self):
        super().setUp()
        self.path = os.path.join(self.data_dir, "enrich2", "dummy.h5")
        self.wt = "ACT"
        self.wt_aa = AA_CODES[CODON_TABLE[self.wt]]
        self.enrich2 = enrich2.Enrich2(self.path, wt_sequence=self.wt)

    def test_parses_nt_variants_into_multi_variant(self):
        nt, _ = self.enrich2.parse_mixed_variant(
            "c.1A>T (p.Thr1Tyr), c.2C>A (p.Thr1Tyr)"
        )
        self.assertEqual(nt, "c.[1A>T;2C>A]")
        #self.assertIsNotNone(hgvsp.multi_variant_re.fullmatch(nt))
        self.assertIsNotNone(re.fullmatch(dna.dna_multi_variant, nt))

    def test_parses_pro_variants_into_multi_variant(self):
        # TODO
        self.enrich2.wt_sequence = "ACTCAA"
        _, pro = self.enrich2.parse_mixed_variant(
            "c.1A>T (p.Thr1Pro), c.4C>A (p.Gln2Lys)"
        )
        self.assertEqual(pro, "p.[Thr1Pro;Gln2Lys]")
        #self.assertIsNotNone(hgvsp.multi_variant_re.fullmatch(pro))
        self.assertIsNotNone(protein.pro_multi_variant.fullmatch(pro))

    def test_error_list_len_different(self):
        # ValueError when attempting tuple unpacking
        with self.assertRaises(ValueError):
            self.enrich2.parse_mixed_variant("c.1A>G, c.2T>A (p.Lys1Arg)")
        with self.assertRaises(ValueError):
            self.enrich2.parse_mixed_variant("p.Lys4Arg, c.2T>A (p.Lys1Arg)")

    def test_variant_order_maintained(self):
        self.enrich2.wt_sequence = "AAAAAT"
        nt, pro = self.enrich2.parse_mixed_variant(
            "c.1= (p.Lys1Ile), c.6T>G (p.Asn2Lys), c.2A>T (p.Lys1Ile)"
        )
        self.assertEqual(nt, "c.[1=;6T>G;2A>T]")
        self.assertEqual(pro, "p.[Lys1Ile;Asn2Lys]")

    @patch.object(
        enrich2.Enrich2, "infer_silent_aa_substitution", return_value="p.Lys1="
    )
    def test_groups_codons(self, patch):
        self.enrich2.wt_sequence = "AAAAAT"
        variant = "c.1= (p.=), c.6T>G (p.Asn2Lys), c.2= (p.=)"
        _, _ = self.enrich2.parse_mixed_variant(variant)
        patch.assert_called_with(*(["c.1=", "c.2="], variant))

    @patch.object(enrich2.Enrich2, "infer_silent_aa_substitution", return_value="p.Lys1=")
    def test_calls_infer_with_synonymous_variants_only(self, patch):
        self.enrich2.wt_sequence = "AAAAAT"
        variant = "c.1= (p.=), c.6T>G (p.Asn2Lys), c.2= (p.Lys1=)"
        _, _ = self.enrich2.parse_mixed_variant(variant)
        patch.assert_called_with(*(["c.1="], variant))

    def test_nt_variant_is_none_special_variant_is_from_synonymous_table(self):
        self.assertEqual(
            (None, constants.enrich2_synonymous),
            self.enrich2.parse_row(
                (constants.enrich2_synonymous, constants.synonymous_table)
            ),
        )

    def test_valueerror_multiple_prefix_types(self):
        with self.assertRaises(ValueError):
            self.enrich2.parse_mixed_variant("c.1A>G (p.=), r.2u>a (p.Lys4Arg)")
        with self.assertRaises(ValueError):
            self.enrich2.parse_mixed_variant("c.1A>G (p.=), n.2T>A (p.Lys4Arg)")
        with self.assertRaises(ValueError):
            self.enrich2.parse_mixed_variant("c.1A>G (p.=), (p.Lys4Arg)")
        with self.assertRaises(ValueError):
            self.enrich2.parse_mixed_variant("c.1A>G (p.=), p.Lys4Arg")
        with self.assertRaises(ValueError):
            self.enrich2.parse_mixed_variant("c.1A>G (p.=), c.2T>A (g.Lys4Arg)")

    def test_doesnt_collapse_single_variants_into_multivariant(self):
        # TODO
        nt, pro = self.enrich2.parse_mixed_variant("c.3T>C (p.=)")
        self.assertEqual(nt, "c.3T>C")
        self.assertEqual(pro, "p.Thr1=")
        #self.assertIsNotNone(hgvsp.single_variant_re.fullmatch(nt))
        self.assertIsNotNone(re.fullmatch(dna.dna_single_variant, nt))
        #self.assertIsNotNone(hgvsp.single_variant_re.fullmatch(pro))
        self.assertIsNotNone(re.fullmatch(protein.pro_multi_variant, pro))

    def test_protein_set_as_nt_when_table_is_not_syn_and_variant_is_special(self):
        nt, pro = self.enrich2.parse_mixed_variant("_wt")
        self.assertEqual(nt, "_wt")
        self.assertEqual(pro, "_wt")

        nt, pro = self.enrich2.parse_mixed_variant("_sy")
        self.assertEqual(nt, "_sy")
        self.assertEqual(pro, "_sy")


class TestInferSilentAASub(ProgramTestCase):
    def setUp(self):
        super().setUp()
        self.path = os.path.join(self.data_dir, "enrich2", "dummy.h5")
        self.enrich2 = enrich2.Enrich2(self.path, wt_sequence="AAA", offset=0)

    def test_valueerror_not_a_sub_event(self):
        with self.assertRaises(exceptions.InvalidVariantType):
            self.enrich2.infer_silent_aa_substitution("c.100_102del")

    def test_index_error_variant_pos_is_out_of_bounds_relative_to_wt(self):
        with self.assertRaises(IndexError):
            self.enrich2.infer_silent_aa_substitution("c.100A>G")

    def test_valueerror_base_in_wt_does_not_match_base_in_hgvs(self):
        self.enrich2.wt_sequence = "TGA"
        with self.assertRaises(ValueError):
            self.enrich2.infer_silent_aa_substitution("c.1A>G")

    def test_correct_wt_aa_inferred(self):
        self.enrich2.wt_sequence = "TCT"
        self.assertEqual("p.Ser1=", self.enrich2.infer_silent_aa_substitution("c.3T>C"))

    def test_correct_aa_position_inferred(self):
        self.enrich2.wt_sequence = "AAAGGGTCT"
        self.assertEqual("p.Ser3=", self.enrich2.infer_silent_aa_substitution("c.9T>C"))

    def test_error_mutant_codon_does_not_match_wild_type(self):
        self.enrich2.wt_sequence = "ATG"
        with self.assertRaises(ValueError):
            self.enrich2.infer_silent_aa_substitution("c.1A>C")

    def test_correctly_infers_aa_from_codon_group(self):
        self.enrich2.wt_sequence = "TTA"
        group = ["c.1T>C", "c.2=", "c.3A>T"]
        self.assertEqual("p.Leu1=", self.enrich2.infer_silent_aa_substitution(group))

    def test_valueerror_mixed_codons_in_group(self):
        with self.assertRaises(ValueError):
            self.enrich2.infer_silent_aa_substitution(["c.1T>C", "c.5T>C"])

    def test_correctly_infers_aa_from_silent_variants(self):
        self.enrich2.wt_sequence = "TTA"
        group = ["c.1=", "c.2=", "c.3="]
        self.assertEqual("p.Leu1=", self.enrich2.infer_silent_aa_substitution(group))


class TestApplyOffset(ProgramTestCase):
    def test_mixed_variant_uses_nt_position_to_compute_codon_pos(self):
        variant = "c.-9A>T (p.Thr2Pro), c.-6C>A (p.Gln3Lys)"
        offset = -10
        self.assertEqual(
            "c.1A>T (p.Thr1Pro), c.4C>A (p.Gln2Lys)",
            enrich2.apply_offset(variant, offset),
        )

    def test_error_position_after_offset_non_positive(self):
        with self.assertRaises(ValueError):
            enrich2.apply_offset("c.1A>T", 10)

        with self.assertRaises(ValueError):
            enrich2.apply_offset("p.Leu1=", 10)

    def test_applies_offset_to_non_mixed_variant(self):
        # TODO
        variant = "n.-455T>A, n.-122A>T, n.-101A>T, n.-42T>A"
        offset = -456
        self.assertEqual(
            "n.1T>A, n.334A>T, n.355A>T, n.414T>A",
            enrich2.apply_offset(variant, offset),
        )
        self.assertEqual("n.1T>A", enrich2.apply_offset("n.-455T>A", offset))

    def test_applies_offset_to_protein_variant_modulo_3(self):
        # TODO
        variant = "p.Leu10=, p.Leu13="
        offset = 10
        self.assertEqual("p.Leu7=, p.Leu10=", enrich2.apply_offset(variant, offset))
        self.assertEqual("p.Leu7=", enrich2.apply_offset("p.Leu10=", offset))

    @patch.object(enrich2.base.BaseProgram, "validate_against_wt_sequence")
    def test_validates_against_wt_sequence(self, patch):
        variant = "c.-9C>T"
        path = os.path.join(self.data_dir, "enrich2", "dummy.h5")
        p = enrich2.Enrich2(path, wt_sequence="ACT")
        enrich2.apply_offset(variant, offset=-10, enrich2=p)  # pass
        patch.assert_called_with(*("c.1C>T",))

    def test_value_error_base_mismatch_after_offset_applied(self):
        variant = "c.-9G>T"
        path = os.path.join(self.data_dir, "enrich2", "dummy.h5")
        p = enrich2.Enrich2(path, wt_sequence="ACT")
        with self.assertRaises(ValueError):
            enrich2.apply_offset(variant, offset=-10, enrich2=p)

    @patch.object(enrich2.base.BaseProgram, "validate_against_protein_sequence")
    def test_validates_against_pro_sequence(self, patch):
        variant = "p.Gly3Leu"
        path = os.path.join(self.data_dir, "enrich2", "dummy.h5")
        p = enrich2.Enrich2(path, wt_sequence="ACG")
        enrich2.apply_offset(variant, offset=6, enrich2=p)  # pass
        patch.assert_called_with(*("p.Gly1Leu",))

    def test_value_error_pro_mismatch_after_offset_applied(self):
        variant = "p.Gly3Leu"
        path = os.path.join(self.data_dir, "enrich2", "dummy.h5")
        p = enrich2.Enrich2(path, wt_sequence="ACG")
        with self.assertRaises(ValueError):
            enrich2.apply_offset(variant, offset=6, enrich2=p)


class TestEnrich2Init(ProgramTestCase):
    def setUp(self):
        super().setUp()
        self.path = os.path.join(self.data_dir, "enrich2", "enrich2.tsv")

    def test_error_is_coding_and_offset_not_mult_of_three(self):
        with self.assertRaises(ValueError):
            enrich2.Enrich2(src=self.path, wt_sequence="ATC", is_coding=True, offset=1)

    def test_ok_is_coding_false_and_offset_not_mult_of_three(self):
        enrich2.Enrich2(src=self.path, wt_sequence="ATC", is_coding=False, offset=1)

    def test_ok_is_coding_and_offset_mult_of_three(self):
        enrich2.Enrich2(src=self.path, wt_sequence="ATC", is_coding=True, offset=-3)


if __name__ == "__main__":
    unittest.main()