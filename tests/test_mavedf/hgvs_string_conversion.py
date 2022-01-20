from unittest import TestCase
from mavetools.mavedf.legacy_to_mave_new import legacy_to_mave_hgvs_nt


class Test(TestCase):
    def test_no_change(self):
        # test no change
        target_seq = "CAATTTGGTTGGTCTGCTAATATGGAA"
        self.assertEqual(legacy_to_mave_hgvs_nt("c.[1=;2=;3=]"), "_wt")
        self.assertEqual(legacy_to_mave_hgvs_nt("n.2="), "_wt")

    def test_one_base_change(self):
        # test one base change
        target_seq = "CAATTTGGTTGGTCTGCTAATATGGAA"
        self.assertEqual(legacy_to_mave_hgvs_nt("c.[1C>A;2=;3=]"), "c.1C>A")
        self.assertEqual(legacy_to_mave_hgvs_nt("c.[4=;5T>A;6=]"), "c.5T>A")
        self.assertEqual(legacy_to_mave_hgvs_nt("c.[4=;5=;6T>A]"), "c.6T>A")
        self.assertEqual(legacy_to_mave_hgvs_nt("c.[4=;5=;6T>G]"), "c.6T>G")

    def test_two_base_change(self):
        # test two base change
        target_seq = "CAATTTGGTTGGTCTGCTAATATGGAA"
        self.assertEqual(legacy_to_mave_hgvs_nt("c.[4T>G;5T>G;6=]"), "c.4_5delinsGG")
        self.assertEqual(legacy_to_mave_hgvs_nt("c.[7G>A;8=;9T>A]"), "c.[7G>A;9T>A]")
        self.assertEqual(legacy_to_mave_hgvs_nt("c.[4T>G;5=;6T>G]"), "c.[4T>G;6T>G]")
        self.assertEqual(legacy_to_mave_hgvs_nt("c.[7=;8G>A;9T>A]"), "c.8_9delinsAA")

    def test_three_base_change(self):
        # test three base change
        target_seq = "CAATTTGGTTGGTCTGCTAATATGGAA"
        self.assertEqual(legacy_to_mave_hgvs_nt("c.[7G>C;8G>T;9T>C]"), "c.7_9delinsCTC")

    def test_four_base_change(self):
        # test four base change
        target_seq = "CAATTTGGTTGGTCTGCTAATATGGAA"
        self.assertEqual(
            legacy_to_mave_hgvs_nt("c.[1C>A;7G>C;8G>T;9T>C]"), "c.[1C>A;7_9delinsCTC]"
        )
        self.assertEqual(
            legacy_to_mave_hgvs_nt("c.[1C>A;2A>C;3A>T;4T>C]"), "c.[1_3delinsACT;4T>C]"
        )
        self.assertEqual(
            legacy_to_mave_hgvs_nt("c.[1C>A;2A>C;8G>T;9T>C]"),
            "c.[1_2delinsAC;8_9delinsTC]",
        )

    def test_five_base_change(self):
        # test five base change
        target_seq = "CAATTTGGTTGGTCTGCTAATATGGAA"
        self.assertEqual(
            legacy_to_mave_hgvs_nt("c.[1C>A;7G>C;8G>T;9T>C;10T>G]"),
            "c.[1C>A;7_9delinsCTC;10T>G]",
        )
