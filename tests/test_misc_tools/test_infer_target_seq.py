from unittest import TestCase
from mavetools.misc_tools import infer_target_seq
from mavehgvs.variant import Variant


class Test(TestCase):
    def test_single_subs(self):
        target_seq = "CAATTTGGTTGGTCTGCTAATATGGAA"
        variant_list = [Variant("n.1C>A"), Variant("n.2A>C"), Variant("n.3A>C"), Variant("n.4T>A"),
                        Variant("n.5T>A"), Variant("n.6T>A"), Variant("n.7G>A"), Variant("n.8G>A"),
                        Variant("n.9T>A"), Variant("n.10T>A"), Variant("n.11G>A"), Variant("n.12G>A"),
                        Variant("n.13T>A"), Variant("n.14C>A"), Variant("n.15T>A"), Variant("n.16G>A"),
                        Variant("n.17C>A"), Variant("n.18T>A"), Variant("n.19A>C"), Variant("n.20A>C"),
                        Variant("n.21T>A"), Variant("n.22A>C"), Variant("n.23T>A"), Variant("n.24G>A"),
                        Variant("n.25G>A"), Variant("n.26A>C"), Variant("n.27A>C")]

        self.assertEqual(infer_target_seq.infer_target_seq(variant_list), target_seq)

    def test_multi_var_single_sub(self):
        target_seq = "CAATTTGGTTGGTCTGCTAATATGGAA"
        variant_list = [Variant("n.[1C>A;3A>C]"), Variant("n.2A>C"), Variant("n.4T>A"),
                        Variant("n.5T>A"), Variant("n.6T>A"), Variant("n.7G>A"), Variant("n.8G>A"),
                        Variant("n.9T>A"), Variant("n.10T>A"), Variant("n.11G>A"), Variant("n.12G>A"),
                        Variant("n.13T>A"), Variant("n.14C>A"), Variant("n.15T>A"), Variant("n.16G>A"),
                        Variant("n.17C>A"), Variant("n.18T>A"), Variant("n.19A>C"), Variant("n.20A>C"),
                        Variant("n.21T>A"), Variant("n.22A>C"), Variant("n.23T>A"), Variant("n.24G>A"),
                        Variant("n.25G>A"), Variant("n.26A>C"), Variant("n.27A>C")]

        self.assertEqual(infer_target_seq.infer_target_seq(variant_list), target_seq)

        variant_list = [Variant("n.1C>A"), Variant("n.2A>C"), Variant("n.3A>C"), Variant("n.4T>A"),
                        Variant("n.5T>A"), Variant("n.6T>A"), Variant("n.7G>A"), Variant("n.8G>A"),
                        Variant("n.9T>A"), Variant("n.10T>A"), Variant("n.11G>A"), Variant("n.12G>A"),
                        Variant("n.13T>A"), Variant("n.14C>A"), Variant("n.15T>A"), Variant("n.16G>A"),
                        Variant("n.17C>A"), Variant("n.18T>A"), Variant("n.19A>C"), Variant("n.20A>C"),
                        Variant("n.21T>A"), Variant("n.22A>C"), Variant("n.23T>A"), Variant("n.24G>A"),
                        Variant("n.26A>C"), Variant("n.[25G>A;27A>C]")]

        self.assertEqual(infer_target_seq.infer_target_seq(variant_list), target_seq)

    def test_incomplete_data(self):
        target_seq = "CAATTTGGTTGGTCTGCTAATATGGAA"
        variant_list = [Variant("n.1C>A"), Variant("n.2A>C"), Variant("n.3A>C"), Variant("n.4T>A"),
                        Variant("n.5T>A"), Variant("n.6T>A"), Variant("n.7G>A"), Variant("n.8G>A"),
                        Variant("n.9T>A"), Variant("n.10T>A"), Variant("n.11G>A"), Variant("n.12G>A"),
                        Variant("n.13T>A"), Variant("n.14C>A"), Variant("n.15T>A"), Variant("n.16G>A"),
                        Variant("n.17C>A"), Variant("n.18T>A"), Variant("n.19A>C"), Variant("n.20A>C"),
                        Variant("n.21T>A"), Variant("n.22A>C"), Variant("n.23T>A"), Variant("n.24G>A"),
                        Variant("n.25G>A"), Variant("n.26A>C")] #, Variant("n.27A>C")]

        self.assertEqual(infer_target_seq.infer_target_seq(variant_list), target_seq)

        target_seq = "CAATTTGGTTGGTCTGCTAATATGGAA"
        variant_list = [Variant("n.1C>A"), Variant("n.2A>C"), Variant("n.3A>C"), Variant("n.4T>A"),
                        #Variant("n.5T>A"), Variant("n.6T>A"), Variant("n.7G>A"), Variant("n.8G>A"),
                        Variant("n.9T>A"), Variant("n.10T>A"), Variant("n.11G>A"), Variant("n.12G>A"),
                        Variant("n.13T>A"), Variant("n.14C>A"), Variant("n.15T>A"), Variant("n.16G>A"),
                        Variant("n.17C>A"), Variant("n.18T>A"), Variant("n.19A>C"), Variant("n.20A>C"),
                        Variant("n.21T>A"), Variant("n.22A>C"), Variant("n.23T>A"), Variant("n.24G>A"),
                        Variant("n.25G>A"), Variant("n.26A>C"), Variant("n.27A>C")]

        self.assertEqual(infer_target_seq.infer_target_seq(variant_list), target_seq)
