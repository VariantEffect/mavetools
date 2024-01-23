import unittest
from mavetools.util.codon import codon_sub_to_variant


class TestCodonSubToMaveHgvs(unittest.TestCase):
    def test_valid_substitution(self):
        valid_cases = [
            (("ATG", "TTG", 1), "c.1A>T"),
            (("ATG", "AAG", 1), "c.2T>A"),
            (("ATG", "ATT", 1), "c.3G>T"),
            (("ATG", "TTG", 2), "c.4A>T"),
            (("ATG", "AAG", 11), "c.32T>A"),
            (("ATG", "ATT", 88), "c.264G>T"),
            (("ATG", "TTT", 1), "c.[1A>T;3G>T]"),
            (("ATG", "TAG", 1, False), "c.[1A>T;2T>A]"),
            (("ATG", "AGT", 1, False), "c.[2T>G;3G>T]"),
            (("ATG", "TAG", 1, True), "c.1_2delinsTA"),
            (("ATG", "AGT", 1, True), "c.2_3delinsGT"),
            (("ATG", "TTT", 2), "c.[4A>T;6G>T]"),
            (("ATG", "TAG", 11, False), "c.[31A>T;32T>A]"),
            (("ATG", "AGT", 88, False), "c.[263T>G;264G>T]"),
            (("ATG", "TAG", 11, True), "c.31_32delinsTA"),
            (("ATG", "AGT", 88, True), "c.263_264delinsGT"),
            (("ATG", "CAT", 1), "c.1_3delinsCAT"),
            (("ATG", "CAT", 101), "c.301_303delinsCAT"),
        ]
        for t, s in valid_cases:
            with self.subTest(t=t, s=s):
                v = codon_sub_to_variant(*t)
                self.assertEqual(s, str(v))

    def test_identical_codons(self):
        valid_cases = [(("ATG", "ATG", 1), "c.="), (("ATG", "ATG", 11), "c.=")]
        for t, s in valid_cases:
            with self.subTest(t=t, s=s):
                v = codon_sub_to_variant(*t)
                self.assertEqual(s, str(v))

    def test_invalid_codon_length(self):
        invalid_cases = [("AT", "TTG", 1), ("ATG", "TTTT", 2), ("ATGG", "CATA", 101)]
        for t in invalid_cases:
            with self.subTest(t=t):
                with self.assertRaises(ValueError):
                    codon_sub_to_variant(*t)

    def test_invalid_codon_bases(self):
        invalid_cases = [("NTG", "TTG", 1), ("ATU", "TTT", 2), ("ATG", "CNT", 101)]
        for t in invalid_cases:
            with self.subTest(t=t):
                with self.assertRaises(ValueError):
                    codon_sub_to_variant(*t)

    def test_invalid_position(self):
        invalid_cases = [
            ("ATG", "TTG", 0),
            ("ATG", "TTT", -12),
            ("ATG", "CAT", "AAA"),
            ("ATG", "AGT", "YFG"),
        ]
        for t in invalid_cases:
            with self.subTest(t=t):
                with self.assertRaises(ValueError):
                    codon_sub_to_variant(*t)


if __name__ == "__main__":
    unittest.main()
