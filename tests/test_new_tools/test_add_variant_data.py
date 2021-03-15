from unittest import TestCase
from mavetools.new_tools.add_variant_data import add_variant_data
from mavetools.new_tools.add_variant_data import parse_additional_hgvs_format


class Test(TestCase):

    # scenario 1 - call function on urn mavedb 00000011-a-1_scores

    # typical case - first argument valid
    def test_add_variant_data(self):
        target_seq = "CAATTTGGTTGGTCTGCTAATATGGAA"
        results = add_variant_data("urn mavedb 00000011-a-1_scores.csv", target_seq)
        # no bases of codon changed
        self.assertEquals(results["variant_codon"][26], "TTT")
        # first base of codon changed
        self.assertEquals(results["variant_codon"][0], "AAA")
        # second base of codon changed
        self.assertEquals(results["variant_codon"][21], "TGT")
        # third base of codon changed
        self.assertEquals(results["variant_codon"][25], "TTG")
        # first and second bases of codon changed
        self.assertEquals(results["variant_codon"][23], "ATA")
        # first and third bases of codon changed
        self.assertEquals(results["variant_codon"][8], "GTG")
        # second and third bases of codon changed
        self.assertEquals(results["variant_codon"][11], "TAC")
        # all three bases of codon changed
        self.assertEquals(results["variant_codon"][12], "AGT")

    # scenario 2 - call function on urn mavedb 00000054-a-1_scores

    # typical case - first argument valid
    def test_add_variant_data2(self):
        target_seq = "ATGACAGCCATCATCAAAGAGATCGTTAGCAGAAACAAAAGGAGATATCAAGAGGATGGATTCGAC" \
                     "TTAGACTTGACCTATATTTATCCAAACATTATTGCTATGGGATTTCCTGCAGAAAGACTTGAAGGC" \
                     "GTATACAGGAACAATATTGATGATGTAGTAAGGTTTTTGGATTCAAAGCATAAAAACCATTACAAG" \
                     "ATATACAATCTTTGTGCTGAAAGACATTATGACACCGCCAAATTTAATTGCAGAGTTGCACAATAT" \
                     "CCTTTTGAAGACCATAACCCACCACAGCTAGAACTTATCAAACCCTTTTGTGAAGATCTTGACCAA" \
                     "TGGCTAAGTGAAGATGACAATCATGTTGCAGCAATTCACTGTAAAGCTGGAAAGGGACGAACTGGT" \
                     "GTAATGATATGTGCATATTTATTACATCGGGGCAAATTTTTAAAGGCACAAGAGGCCCTAGATTTC" \
                     "TATGGGGAAGTAAGGACCAGAGACAAAAAGGGAGTAACTATTCCCAGTCAGAGGCGCTATGTGTAT" \
                     "TATTATAGCTACCTGTTAAAGAATCATCTGGATTATAGACCAGTGGCACTGTTGTTTCACAAGATG" \
                     "ATGTTTGAAACTATTCCAATGTTCAGTGGCGGAACTTGCAATCCTCAGTTTGTGGTCTGCCAGCTA" \
                     "AAGGTGAAGATATATTCCTCCAATTCAGGACCCACACGACGGGAAGACAAGTTCATGTACTTTGAG" \
                     "TTCCCTCAGCCGTTACCTGTGTGTGGTGATATCAAAGTAGAGTTCTTCCACAAACAGAACAAGATG" \
                     "CTAAAAAAGGACAAAATGTTTCACTTTTGGGTAAATACATTCTTCATACCAGGACCAGAGGAAACC" \
                     "TCAGAAAAAGTAGAAAATGGAAGTCTATGTGATCAAGAAATCGATAGCATTTGCAGTATAGAGCGT" \
                     "GCAGATAATGACAAGGAATATCTAGTACTTACTTTAACAAAAAATGATCTTGACAAAGCAAATAAA" \
                     "GACAAAGCCAACCGATACTTTTCTCCAAATTTTAAGGTGAAGCTGTACTTCACAAAAACAGTAGAG" \
                     "GAGCCGTCAAATCCAGAGGCTAGCAGTTCAACTTCTGTAACACCAGATGTTAGTGACAATGAACCT" \
                     "GATCATTATAGATATTCTGACACCACTGACTCTGATCCAGAGAATGAACCTTTTGATGAAGATCAG" \
                     "CATACACAAATTACAAAAGTCTGA"
        results = add_variant_data("urn mavedb 00000054-a-1_scores.csv", target_seq)
        # no bases of codon changed
        self.assertEquals(results["variant_codon"][8676], None)
        # first base of codon changed
        self.assertEquals(results["variant_codon"][8677], "TTG")
        # second base of codon changed
        self.assertEquals(results["variant_codon"][8643], "AGA")
        # third base of codon changed
        self.assertEquals(results["variant_codon"][8678], "ACT")
        # first and second bases of codon changed
        self.assertEquals(results["variant_codon"][8673], "TGG")
        # first and third bases of codon changed
        self.assertEquals(results["variant_codon"][8675], "TTT")
        # second and third bases of codon changed
        self.assertEquals(results["variant_codon"][8671], "ACT")
        # all three bases of codon changed
        self.assertEquals(results["variant_codon"][8674], "TAT")
