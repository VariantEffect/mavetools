from unittest import TestCase
from mavetools.misc_tools import delins_to_subs


class Test(TestCase):
    def test_two_base_delins(self):
        target_seq = "CAATTTGGTTGGTCTGCTAATATGGAA"
        self.assertEqual(delins_to_subs.delins_to_subs(target_seq, "c.4_5delinsGG"), "c.[4T>G;5T>G]")

    def test_base_identical_to_target(self):
        target_seq = "CAATTTGGTTGGTCTGCTAATATGGAA"
        self.assertEqual(delins_to_subs.delins_to_subs(target_seq, "c.4_6delinsGTG"), "c.[4T>G;6T>G]")

    def test_offset(self):
        target_seq = "CAATTTGGTTGGTCTGCTAATATGGAA"
        self.assertEqual(delins_to_subs.delins_to_subs(target_seq, "c.1_3delinsGTG", offset=3), "c.[1T>G;3T>G]")
        self.assertEqual(delins_to_subs.delins_to_subs(target_seq, "c.1_2delinsGT", offset=25), "c.[1A>G;2A>T]")

    def test_variant_identical_to_target(self):
        target_seq = "CAATTTGGTTGGTCTGCTAATATGGAA"
        with self.assertRaises(ValueError):
            delins_to_subs.delins_to_subs(target_seq, "c.1_3delinsGTG", offset=25)

    def test_variant_out_of_bounds(self):
        target_seq = "CAATTTGGTTGGTCTGCTAATATGGAA"
        with self.assertRaises(ValueError):
            delins_to_subs.delins_to_subs(target_seq, "c.1_3delinsTTT", offset=3)
