import os
import unittest
from unittest.mock import patch

from mavehgvs.exceptions import MaveHgvsParseError
from mavetools.mavedbconvert import base, exceptions

from tests import ProgramTestCase


class BaseTest(base.BaseProgram):
    """
    Class with no I/O or parsing functionality for testing purposes.
    Required because base.BaseProgram is an abstract base class.
    """

    def __init__(
        self,
        src,
        wt_sequence,
        offset=0,
        dst=None,
        one_based=False,
        skip_header_rows=0,
        skip_footer_rows=0,
        score_column=None,
        hgvs_column=None,
        input_type=None,
        sheet_name=None,
        is_coding=True,
    ):
        super().__init__(
            src=src,
            wt_sequence=wt_sequence,
            offset=offset,
            dst=dst,
            one_based=one_based,
            skip_header_rows=skip_header_rows,
            skip_footer_rows=skip_footer_rows,
            sheet_name=sheet_name,
            score_column=score_column,
            input_type=input_type,
            hgvs_column=hgvs_column,
            is_coding=is_coding,
        )

    def load_input_file(self):
        pass

    def parse_row(self, row):
        pass

    def parse_input(self, df):
        pass


class TestPaths(ProgramTestCase):
    """
    Test __init__ correctly sets up read and write directories,
    etc.
    """

    def setUp(self):
        super().setUp()
        self.src = os.path.join(self.data_dir, "enrich", "enrich.tsv")
        self.src_with_spaces = os.path.join(self.data_dir, "enrich", "enrich   .tsv")
        self.h5_src = os.path.join(self.data_dir, "enrich2", "dummy.h5")

    def test_sets_directory_as_input_directory_if_dst_is_none(self):
        p = BaseTest(src=self.src, dst=None, wt_sequence="AAA")
        self.assertEqual(p.dst, os.path.join(self.data_dir, "enrich"))

    def test_error_file_not_readable(self):
        with self.assertRaises(IOError):
            BaseTest(src="a file", dst=None, wt_sequence="AAA")

    def test_expands_user_and_norms_dst(self):
        p = BaseTest(src=self.src, dst="~//user//", wt_sequence="AAA")
        self.assertEqual(p.dst, os.path.join(os.path.expanduser("~"), "user"))

    def test_dir_with_input_fname_appended_when_h5_and_dst_is_none(self):
        p = BaseTest(src=self.h5_src, dst=None, wt_sequence="AAA")
        self.assertEqual(p.dst, os.path.join(self.data_dir, "enrich2", "dummy"))

    def test_creates_directory_tree_if_it_doesnt_exist(self):
        output = os.path.join(self.data_dir, "enrich2", "outer_dir", "inner_dir")
        BaseTest(src=self.h5_src, dst=output, wt_sequence="AAA")
        self.assertTrue(os.path.isdir(output))

    @patch("os.access")
    def test_checks_read_permission(self, patch):
        p = BaseTest(src=self.src, dst=None, wt_sequence="AAA")
        self.assertEqual(patch.call_args_list[0][0], (p.src, os.R_OK))

    @patch("os.access")
    def test_checks_write_permission(self, patch):
        p = BaseTest(src=self.src, dst=None, wt_sequence="AAA")
        self.assertEqual(patch.call_args_list[1][0], (p.dst, os.W_OK))

    def test_splits_src_into_filename_and_ext(self):
        p = BaseTest(src=self.src, dst=None, wt_sequence="AAA")
        self.assertEqual(p.src_filename, "enrich")
        self.assertEqual(p.ext, ".tsv")

    def test_lower_cases_ext(self):
        p = BaseTest(src=self.src.replace("tsv", "TSV"), wt_sequence="AAA")
        self.assertEqual(p.ext, ".tsv")

    def test_dst_filename_replaces_whitespace_with_underscores(self):
        p = BaseTest(src=self.src_with_spaces, wt_sequence="AAA")
        self.assertEqual(p.dst_filename, "mavedb_enrich_.csv")

    def test_output_file_joins_dst_and_dst_filename(self):
        p = BaseTest(src=self.src, wt_sequence="AAA")
        self.assertEqual(
            p.output_file, os.path.join(self.data_dir, "enrich", "mavedb_enrich.csv")
        )

    def test_output_directory_expands_user_and_norms_path(self):
        p = BaseTest(src=self.src, wt_sequence="AAA")
        p.dst = "~//user//"
        self.assertEqual(
            p.output_directory, os.path.join(os.path.expanduser("~"), "user")
        )


class TestWtSequence(ProgramTestCase):
    """
    Test __init__ correctly sets up sequence information etc.
    """

    def setUp(self):
        super().setUp()
        self.src = os.path.join(self.data_dir, "enrich", "enrich.tsv")
        self.src_with_spaces = os.path.join(self.data_dir, "enrich", "enrich   .tsv")
        self.h5_src = os.path.join(self.data_dir, "enrich2", "dummy.h5")

    def test_value_error_coding_offset_not_multiple_of_three(self):
        with self.assertRaises(ValueError):
            BaseTest(src=self.src, wt_sequence="ATCA", offset=-1)

    # --- Test property setters --- #
    def test_wt_setter_upper_cases_wt_sequence(self):
        p = BaseTest(src=self.src, wt_sequence="AAA")
        p.wt_sequence = "ggg"
        self.assertEqual(p.wt_sequence, "GGG")

    def test_wt_setter_uses_full_wt_sequence_and_ignores_offset(self):
        p = BaseTest(src=self.src, wt_sequence="AAA")
        p.offset = 3
        p.wt_sequence = "ATGCGA"
        self.assertEqual(p.wt_sequence, "ATGCGA")

    def test_wt_setter_creates_codons_from_wt_sequence_and_ignores_offset(self):
        p = BaseTest(src=self.src, wt_sequence="AAA")
        p.offset = 3
        p.wt_sequence = "ATGCGA"
        self.assertEqual(p.codons, ["ATG", "CGA"])

    def test_wt_setter_value_error_not_valid_wt_sequence(self):
        p = BaseTest(src=self.src, wt_sequence="AAA")
        with self.assertRaises(ValueError):
            p.wt_sequence = "fff"

    # def test_error_is_coding_and_offset_not_multiple_of_three(self):
    #     with self.assertRaises(ValueError):
    #         BaseTest(
    #             src=self.src,
    #             wt_sequence='AAA',
    #             is_coding=True,
    #             offset=2
    #         )


class TestBaseProgramValidateAgainstWTSeq(ProgramTestCase):
    def setUp(self):
        super().setUp()
        self.src = os.path.join(self.data_dir, "enrich", "enrich.tsv")
        self.base = BaseTest(src=self.src, wt_sequence="ATG", one_based=True)

    def test_error_not_a_dna_sub(self):
        with self.assertRaises(exceptions.InvalidVariantType):
            self.base.validate_against_wt_sequence("p.Gly1Leu")

    def test_passes_on_special_and_silent(self):
        self.base.validate_against_wt_sequence("c.1=")

    def test_throws_error_with_invalid_variants(self):
        with self.assertRaises(MaveHgvsParseError):
            self.base.validate_against_wt_sequence("_wt")
        with self.assertRaises(MaveHgvsParseError):
            self.base.validate_against_wt_sequence("_sy")

    def test_passes_when_reference_base_matches(self):
        self.base.one_based = False
        self.base.validate_against_wt_sequence("c.1T>G")

        self.base.one_based = True
        self.base.validate_against_wt_sequence("c.1A>G")

    def test_error_when_reference_base_doesnt_match_wt(self):
        with self.assertRaises(ValueError):
            self.base.one_based = False
            self.base.validate_against_wt_sequence("c.1C>G")

        with self.assertRaises(ValueError):
            self.base.one_based = True
            self.base.validate_against_wt_sequence("c.1T>G")

    def test_error_negative_position(self):
        with self.assertRaises(MaveHgvsParseError):
            self.base.one_based = True
            self.base.validate_against_wt_sequence("c.0T>G")

    def test_validates_multi(self):
        self.base.validate_against_wt_sequence("c.[1A>G;2T>G]")
        with self.assertRaises(ValueError):
            self.base.validate_against_wt_sequence("c.[1A>G;2A>G]")

    def test_index_error_index_extends_beyond_indexable_wt_seq(self):
        with self.assertRaises(IndexError):
            self.base.one_based = True
            self.base.validate_against_wt_sequence("c.4G>A")

        with self.assertRaises(IndexError):
            self.base.one_based = False
            self.base.validate_against_wt_sequence("c.3G>A")


class TestBaseProgramValidateAgainstProteinSeq(ProgramTestCase):
    def setUp(self):
        super().setUp()
        self.src = os.path.join(self.data_dir, "enrich", "enrich.tsv")
        self.base = BaseTest(src=self.src, wt_sequence="ATGAAA", one_based=True)

    def test_error_not_a_protein_sub(self):
        with self.assertRaises(exceptions.InvalidVariantType):
            self.base.validate_against_protein_sequence("c.1A>G")

    def test_passes_on_special_and_silent(self):
        self.base.validate_against_protein_sequence("p.=")

    def test_passes_when_reference_aa_matches(self):
        self.base.validate_against_protein_sequence("p.Met1Lys")

    def test_error_when_reference_base_doesnt_match_wt(self):
        with self.assertRaises(ValueError):
            self.base.validate_against_protein_sequence("p.Met2Lys")

    def test_error_invalid_position(self):
        # TODO
        with self.assertRaises(ValueError):
            self.base.validate_against_protein_sequence("p.Met0Lys")

    def test_validates_multi(self):
        # TODO
        self.base.validate_against_protein_sequence("p.[Met1Lys;Lys2=]")
        with self.assertRaises(ValueError):
            self.base.validate_against_protein_sequence("p.[Met1Lys;Met2=]")

    def test_error_invalid_position_multi(self):
        # TODO
        with self.assertRaises(ValueError):
            self.base.validate_against_protein_sequence("p.[Met0Lys;Lys2=]")

    def test_index_error_index_extends_beyond_indexable_pro_seq(self):
        self.base.validate_against_protein_sequence("p.Lys2Met")

        with self.assertRaises(IndexError):
            self.base.validate_against_protein_sequence("p.Met3Lys")


if __name__ == "__main__":
    unittest.main()