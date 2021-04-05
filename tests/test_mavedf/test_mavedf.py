from unittest import TestCase
from io import StringIO
from mavetools.mavedf.mavedf import MaveDf


class Test(TestCase):

    # typical case - no bases of codon changed
    def test_no_change(self):
        data = StringIO("# Accession:\n"
                        "# Downloaded (UTC):\n"
                        "# Licence:\n"
                        "# Licence URL:\n"
                        "accession,hgvs_nt,hgvs_pro,score\n"
                        "urn:mavedb:00000011-a-1#27,c.[4=;5=;6=],p.Phe2=,0.0")
        target_seq = "CAATTTGGTTGGTCTGCTAATATGGAA"
        df = MaveDf(data)
        df.add_variant_data(target_seq)
        self.assertEqual(df.pandas_df["variant_codon"][0], None)  # TTT but no way to know after converting from legacy

        data = StringIO("# Accession:\n"
                        "# Downloaded (UTC):\n"
                        "# Licence:\n"
                        "# Licence URL:\n"
                        "accession,hgvs_nt,hgvs_pro,score\n"
                        "urn:mavedb:00000054-a-1#8677,_wt,_wt,1.0")
        target_seq = "ATGACA"
        df = MaveDf(data)
        df.add_variant_data(target_seq)
        self.assertEqual(df.pandas_df["variant_codon"][0], None)

    # typical case - first base of codon changed
    def test_first_base_change(self):
        data = StringIO("# Accession:\n"
                        "# Downloaded (UTC):\n"
                        "# Licence:\n"
                        "# Licence URL:\n"
                        "accession,hgvs_nt,hgvs_pro,score\n"
                        "urn:mavedb:00000011-a-1#1,c.[1C>A;2=;3=],p.Gln1Lys,0.016527088")
        target_seq = "CAATTTGGTTGGTCTGCTAATATGGAA"
        df = MaveDf(data)
        df.add_variant_data(target_seq)
        self.assertEqual(df.pandas_df["variant_codon"][0], "AAA")

        data = StringIO("# Accession:\n"
                        "# Downloaded (UTC):\n"
                        "# Licence:\n"
                        "# Licence URL:\n"
                        "accession,hgvs_nt,hgvs_pro,score\n"
                        "urn:mavedb:00000054-a-1#8678,c.1A>T,p.Met1Leu,-1.545613663")
        target_seq = "ATGACA"
        df = MaveDf(data)
        df.add_variant_data(target_seq)
        self.assertEqual(df.pandas_df["variant_codon"][0], "TTG")

    # typical case - second base of codon changed
    def test_second_base_change(self):
        data = StringIO("# Accession:\n"
                        "# Downloaded (UTC):\n"
                        "# Licence:\n"
                        "# Licence URL:\n"
                        "accession,hgvs_nt,hgvs_pro,score\n"
                        "urn:mavedb:00000011-a-1#22,c.[4=;5T>G;6=],p.Phe2Cys,-0.164530439\n")
        target_seq = "CAATTTGGTTGGTCTGCTAATATGGAA"
        df = MaveDf(data)
        df.add_variant_data(target_seq)
        self.assertEqual(df.pandas_df["variant_codon"][0], "TGT")

        data = StringIO("# Accession:\n"
                        "# Downloaded (UTC):\n"
                        "# Licence:\n"
                        "# Licence URL:\n"
                        "accession,hgvs_nt,hgvs_pro,score\n"
                        "urn:mavedb:00000054-a-1#8644,c.5C>G,p.Thr2Arg,0.551584254")
        target_seq = "ATGACA"
        df = MaveDf(data)
        df.add_variant_data(target_seq)
        self.assertEqual(df.pandas_df["variant_codon"][0], "AGA")

    # typical case - third base of codon changed
    def test_third_base_change(self):
        data = StringIO("# Accession:\n"
                        "# Downloaded (UTC):\n"
                        "# Licence:\n"
                        "# Licence URL:\n"
                        "accession,hgvs_nt,hgvs_pro,score\n"
                        "urn:mavedb:00000011-a-1#26,c.[4=;5=;6T>G],p.Phe2Leu,-0.093330569")
        target_seq = "CAATTTGGTTGGTCTGCTAATATGGAA"
        df = MaveDf(data)
        df.add_variant_data(target_seq)
        self.assertEqual(df.pandas_df["variant_codon"][0], "TTG")

        data = StringIO("# Accession:\n"
                        "# Downloaded (UTC):\n"
                        "# Licence:\n"
                        "# Licence URL:\n"
                        "accession,hgvs_nt,hgvs_pro,score\n"
                        "urn:mavedb:00000054-a-1#8679,c.3G>T,p.Met1Ile,-0.679812995")
        target_seq = "ATGACA"
        df = MaveDf(data)
        df.add_variant_data(target_seq)
        self.assertEqual(df.pandas_df["variant_codon"][0], "ACT")

    # typical case - first and second bases of codon changed
    def test_first_second_base_change(self):
        data = StringIO("# Accession:\n"
                        "# Downloaded (UTC):\n"
                        "# Licence:\n"
                        "# Licence URL:\n"
                        "accession,hgvs_nt,hgvs_pro,score\n"
                        "urn:mavedb:00000011-a-1#24,c.[1C>A;2A>T;3=],p.Gln1Ile,0.016768928")
        target_seq = "CAATTTGGTTGGTCTGCTAATATGGAA"
        df = MaveDf(data)
        df.add_variant_data(target_seq)
        self.assertEqual(df.pandas_df["variant_codon"][0], "ATA")

        data = StringIO("# Accession:\n"
                        "# Downloaded (UTC):\n"
                        "# Licence:\n"
                        "# Licence URL:\n"
                        "accession,hgvs_nt,hgvs_pro,score\n"
                        "urn:mavedb:00000054-a-1#8674,c.1_2delinsTG,p.Met1Trp,-1.300059776")
        target_seq = "ATGACA"
        df = MaveDf(data)
        df.add_variant_data(target_seq)
        self.assertEqual(df.pandas_df["variant_codon"][0], "TGG")

    # typical case - first and third bases of codon changed
    def test_first_third_base_change(self):
        data = StringIO("# Accession:\n"
                        "# Downloaded (UTC):\n"
                        "# Licence:\n"
                        "# Licence URL:\n"
                        "accession,hgvs_nt,hgvs_pro,score\n"
                        "urn:mavedb:00000011-a-1#9,c.[4T>G;5=;6T>G],p.Phe2Val,-0.66351736")
        target_seq = "CAATTTGGTTGGTCTGCTAATATGGAA"
        df = MaveDf(data)
        df.add_variant_data(target_seq)
        self.assertEqual(df.pandas_df["variant_codon"][0], "GTG")

        data = StringIO("# Accession:\n"
                        "# Downloaded (UTC):\n"
                        "# Licence:\n"
                        "# Licence URL:\n"
                        "accession,hgvs_nt,hgvs_pro,score\n"
                        "urn:mavedb:00000054-a-1#8676,c.[1A>T;3G>T],p.Met1Phe,NA")
        target_seq = "ATGACA"
        df = MaveDf(data)
        df.add_variant_data(target_seq)
        self.assertEqual(df.pandas_df["variant_codon"][0], "TTT")

    # typical case - second and third bases of codon changed
    def test_second_third_base_change(self):
        data = StringIO("# Accession:\n"
                        "# Downloaded (UTC):\n"
                        "# Licence:\n"
                        "# Licence URL:\n"
                        "accession,hgvs_nt,hgvs_pro,score\n"
                        "urn:mavedb:00000011-a-1#12,c.[4=;5T>A;6T>C],p.Phe2Tyr,-0.000442088")
        target_seq = "CAATTTGGTTGGTCTGCTAATATGGAA"
        df = MaveDf(data)
        df.add_variant_data(target_seq)
        self.assertEqual(df.pandas_df["variant_codon"][0], "TAC")

        data = StringIO("# Accession:\n"
                        "# Downloaded (UTC):\n"
                        "# Licence:\n"
                        "# Licence URL:\n"
                        "accession,hgvs_nt,hgvs_pro,score\n"
                        "urn:mavedb:00000054-a-1#8672,c.2_3delinsCT,p.Met1Thr,NA")
        target_seq = "ATGACA"
        df = MaveDf(data)
        df.add_variant_data(target_seq)
        self.assertEqual(df.pandas_df["variant_codon"][0], "ACT")

    # typical case - all three bases of codon changed
    def test_all_base_change(self):
        data = StringIO("# Accession:\n"
                        "# Downloaded (UTC):\n"
                        "# Licence:\n"
                        "# Licence URL:\n"
                        "accession,hgvs_nt,hgvs_pro,score\n"
                        "urn:mavedb:00000011-a-1#13,c.[1C>A;2A>G;3A>T],p.Gln1Ser,-0.019056815")
        target_seq = "CAATTTGGTTGGTCTGCTAATATGGAA"
        df = MaveDf(data)
        df.add_variant_data(target_seq)
        self.assertEqual(df.pandas_df["variant_codon"][0], "AGT")

        data = StringIO("# Accession:\n"
                        "# Downloaded (UTC):\n"
                        "# Licence:\n"
                        "# Licence URL:\n"
                        "accession,hgvs_nt,hgvs_pro,score\n"
                        "urn:mavedb:00000054-a-1#8675,c.1_3delinsTAT,p.Met1Tyr,-0.872033698")
        target_seq = "ATGACA"
        df = MaveDf(data)
        df.add_variant_data(target_seq)
        self.assertEqual(df.pandas_df["variant_codon"][0], "TAT")

    # scenario 1 - call function on urn mavedb 00000011-a-1_scores
""""
    # typical case - test drop_accession and ret_meta
    def test_add_variant_data_drop_accession(self):
        target_seq = "CAATTTGGTTGGTCTGCTAATATGGAA"
        results = add_variant_data("urn mavedb 00000011-a-1_scores.csv", target_seq, True, True)
        self.assertEquals(results[0].iat[0, 0], "c.[1C>A;2=;3=]")

    # edge case - invalid second or third argument boolean
    def test_add_variant_data_not_boolean_arg(self):
        with self.assertRaises(TypeError):
            target_seq = "CAATTTGGTTGGTCTGCTAATATGGAA"
            add_variant_data("urn mavedb 00000011-a-1_scores.csv", target_seq, drop_accession="cat")
        with self.assertRaises(TypeError):
            target_seq = "CAATTTGGTTGGTCTGCTAATATGGAA"
            add_variant_data("urn mavedb 00000011-a-1_scores.csv", target_seq, ret_meta="cat")

    # test helper functions

    # test parse_additional_hgvs_format
    def test_add_variant_data_parser(self):
        sub_one, sub_two, sub_three, sub_one_nuc, sub_two_nuc, sub_three_nuc = parse_additional_hgvs_format(
            "c.[1C>A;2=;3=]")
        self.assertEquals(sub_one, 0)
        self.assertEquals(sub_two, None)
        self.assertEquals(sub_three, None)
        self.assertEquals(sub_one_nuc, "A")
        self.assertEquals(sub_two_nuc, None)
        self.assertEquals(sub_three_nuc, None) """