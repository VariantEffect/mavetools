from unittest import TestCase
from mavetools.misc_tools import infer_target_seq
from mavehgvs.variant import Variant


class Test(TestCase):
    def test_target_seq_inference(self):
        target_seq = "CAATTTGGTTGGTCTGCTAATATGGAA"
        variant_list = [Variant("1C>A"), Variant("2A>C"), Variant("3A>C"), Variant("4T>A"), Variant("5T>A"),
                        Variant("6T>A"), Variant("7G>A"), Variant("8G>A"), Variant("9T>A"), Variant("10T>A"),
                        Variant("11G>A"), Variant("12G>A"), Variant("13T>A"), Variant("14C>A"), Variant("15T>A"),
                        Variant("16G>A"), Variant("17C>A"), Variant("18T>A"), Variant("19A>C"), Variant("20A>C"),
                        Variant("21T>A"), Variant("22A>C"), Variant("23T>A"), Variant("24G>A"), Variant("25G>A"),
                        Variant("26A>C"), Variant("27A>C")]

        self.assertEqual(infer_target_seq.infer_target_seq(variant_list), target_seq)
