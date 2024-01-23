import unittest
from mavehgvs import Variant
from mavetools.convert.enrich import seqid_to_variant


class TestSeqidToVariant(unittest.TestCase):
    def test_protein_wt(self):
        wtseq = "MGLYTAKLVEN*"
        seqid_tuples = [
            ("p.Met1Leu", "0-L"),
            ("p.Tyr4Asp", "3-D"),
            ("p.[Met1Leu;Tyr4Asp]", "0,3-L,D"),
            ("p.[Met1Leu;Tyr4Asp;Ter12Asn]", "0,3,11-L,D,N"),
            ("p.Tyr4Ter", "3-*"),
            ("p.Ter12Gly", "11-G"),
        ]

        for hgvs, seqid in seqid_tuples:
            with self.subTest(hgvs=hgvs, seqid=seqid, wtseq=wtseq):
                self.assertEqual(Variant(hgvs), seqid_to_variant(seqid, wtseq))

    def test_invalid_seqids(self):
        invalid_seqids = [
            "MGLYTAKLVEN*",
            "0,3-L",
            "-L",
            "-1-D",
            "3-L,D",
            "0,3;L,D",
            "0;3-L;D",
            "0,3-L,B",
            "0-X",
        ]

        for seqid in invalid_seqids:
            with self.subTest(seqid=seqid):
                with self.assertRaises(ValueError):
                    seqid_to_variant(seqid, "AAA")

    def test_wt_length(self):
        wtseq = "WT"
        valid_seqids = ["3-D", "0,3-L,D", "0,3,12-L,D,N", "3-*"]

        for seqid in valid_seqids:
            with self.subTest(seqid=seqid, wtseq=wtseq):
                with self.assertRaises(ValueError):
                    seqid_to_variant(seqid, wtseq)

    def test_synonymous_not_allowed(self):
        wtseq = "WT*"
        valid_seqids = ["0-W", "0,3-W,D", "0,1-L,T", "2-*"]

        for seqid in valid_seqids:
            with self.subTest(seqid=seqid, wtseq=wtseq):
                with self.assertRaises(ValueError):
                    seqid_to_variant(seqid, wtseq)


if __name__ == "__main__":
    unittest.main()
