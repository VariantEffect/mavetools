from unittest import TestCase
from mavetools.legacy_to_mave_hgvs.legacy_to_mave_hgvs import legacy_to_mave_hgvs_nt

class Test(TestCase):
    def test_legacy_to_mave_hgvs_nt(self):
        target_seq = "CAATTTGGTTGGTCTGCTAATATGGAA"
        # test no change
        self.assertEqual(legacy_to_mave_hgvs_nt("c.[1=;2=;3=]", target_seq), "_wt")
        # test one base change
        self.assertEqual(legacy_to_mave_hgvs_nt("c.[1C>A;2=;3=]", target_seq), "c.1C>A")
        self.assertEqual(legacy_to_mave_hgvs_nt("c.[4=;5T>A;6=]", target_seq), "c.5T>A")
        self.assertEqual(legacy_to_mave_hgvs_nt("c.[4=;5=;6T>A]", target_seq), "c.6T>A")
        self.assertEqual(legacy_to_mave_hgvs_nt("c.[4=;5=;6T>G]", target_seq), "c.6T>G")
        # test two base change
        self.assertEqual(legacy_to_mave_hgvs_nt("c.[4T>G;5T>G;6=]", target_seq), "c.4_5delinsGG")
        self.assertEqual(legacy_to_mave_hgvs_nt("c.[7G>A;8=;9T>A]", target_seq), "c.[7G>A;9T>A]")
        self.assertEqual(legacy_to_mave_hgvs_nt("c.[4T>G;5=;6T>G]", target_seq), "c.[4T>G;6T>G]")
        self.assertEqual(legacy_to_mave_hgvs_nt("c.[7=;8G>A;9T>A]", target_seq), "c.8_9delinsAA")
        # test three base change
        self.assertEqual(legacy_to_mave_hgvs_nt("c.[7G>C;8G>T;9T>C]", target_seq), "c.7_9delinsCTC")

        # self.fail()
