from unittest import TestCase
from mavetools.mavevariant.mavevariant import MaveVariant


class Test(TestCase):
    def test_no_change(self):
        # test no change
        target_seq = "CAATTTGGTTGGTCTGCTAATATGGAA"
        legacy_hgvs = "c.[1=;2=;3=]"
        variant = MaveVariant(legacy_hgvs, target_seq)
        self.assertEqual(variant.mave_hgvs, "_wt")

    def test_one_base_change(self):
        # test one base change
        target_seq = "CAATTTGGTTGGTCTGCTAATATGGAA"
        variant = MaveVariant("c.[1C>A;2=;3=]", target_seq)
        self.assertEqual(variant.mave_hgvs, "c.1C>A")
        variant = MaveVariant("c.[4=;5T>A;6=]", target_seq)
        self.assertEqual(variant.mave_hgvs, "c.5T>A")
        variant = MaveVariant("c.[4=;5=;6T>A]", target_seq)
        self.assertEqual(variant.mave_hgvs, "c.6T>A")
        variant = MaveVariant("c.[4=;5=;6T>G]", target_seq)
        self.assertEqual(variant.mave_hgvs, "c.6T>G")

    def test_two_base_change(self):
        # test two base change
        target_seq = "CAATTTGGTTGGTCTGCTAATATGGAA"
        variant = MaveVariant("c.[4T>G;5T>G;6=]", target_seq)
        self.assertEqual(variant.mave_hgvs, "c.4_5delinsGG")
        variant = MaveVariant("c.[7G>A;8=;9T>A]", target_seq)
        self.assertEqual(variant.mave_hgvs, "c.[7G>A;9T>A]")
        variant = MaveVariant("c.[4T>G;5=;6T>G]", target_seq)
        self.assertEqual(variant.mave_hgvs, "c.[4T>G;6T>G]")
        variant = MaveVariant("c.[7=;8G>A;9T>A]", target_seq)
        self.assertEqual(variant.mave_hgvs, "c.8_9delinsAA")

    def test_three_base_change(self):
        # test three base change
        target_seq = "CAATTTGGTTGGTCTGCTAATATGGAA"
        variant = MaveVariant("c.[7G>C;8G>T;9T>C]", target_seq)
        self.assertEqual(variant.mave_hgvs, "c.7_9delinsCTC")

    # def test_
