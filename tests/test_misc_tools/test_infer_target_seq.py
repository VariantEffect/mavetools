from unittest import TestCase
from mavetools.misc_tools import infer_target_seq


class Test(TestCase):
    def test_target_seq_inference(self):
        target_seq = "CAATTTGGTTGGTCTGCTAATATGGAA"
        variant_list = []
        self.assertEqual(infer_target_seq.infer_target_seq(variant_list), target_seq)
