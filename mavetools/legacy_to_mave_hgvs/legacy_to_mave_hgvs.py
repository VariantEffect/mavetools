from collections import Counter
import re


def legacy_to_mave_hgvs_nt(hgvs_legacy, target_seq):
    """
    This function converts a legacy hgvs_nt formatted string (i.e., c.[1C>A;2=;3=]) and converts it
    to the standard format (i.e., c.1delinsA). If string is already in standard format, the input
    string is returned as-is.

    Parameters
    ----------
    hgvs_legacy (string): legacy format hgvs string
    target_seq (string): target sequence

    Returns
    -------
    standard_hgvs (string): standard format hgvs string
    """
    # determine which bases had a change
    changes = [letter for letter in hgvs_legacy if letter in "=>"]

    # if changes < 3, already in standard_hgvs format
    if len(changes) < 3:
        return hgvs_legacy

    # count instances of changes
    count = Counter(changes)
    count = count[">"]

    if count == 0:
        # variant_codon is wild-type
        standard_hgvs = "_wt"

    elif count == 1:
        # variant_codon has one nucleotide substitution
        # get locations of nucleotide substitutions
        sub = re.split("\[", hgvs_legacy)[1]
        sub = re.split("[a-zA-Z.>;=]", sub)
        # get nucleotides of substitutions
        sub_nuc = hgvs_legacy.split("]")[0]
        sub_nuc = re.split("\[", sub_nuc)[1]
        sub_nuc = re.split("[a-z0-9.>;=]", sub_nuc)
        sub_nuc = [letter for letter in sub_nuc if letter in 'ACTG' and letter != ""]
        if changes[0] == ">" and changes[1] == "=" and changes[2] == "=":
            sub_one = int(sub[0])
        if changes[0] == "=" and changes[1] == ">" and changes[2] == "=":
            sub_one = int(sub[2])
        if changes[0] == "=" and changes[1] == "=" and changes[2] == ">":
            sub_one = int(sub[4])
        # get substituted nucleotide
        sub_one_nuc = sub_nuc[1]
        # construct hgvs_nt string
        standard_hgvs = "c." + str(sub_one) + target_seq[sub_one - 1] + ">" + sub_one_nuc

    elif count == 2:
        # keep track of whether substitutions are adjacent
        adjacent = True
        # variant_codon has two nucleotide substitutions
        # get indices of nucleotide substitutions
        sub = re.split("\[", hgvs_legacy)[1]
        sub = re.split("[a-zA-Z.>;=]", sub)
        # get nucleotides of substitutions
        sub_nuc = hgvs_legacy.split("]")[0]
        sub_nuc = re.split("\[", sub_nuc)[1]
        sub_nuc = re.split("[a-z0-9.>;=]", sub_nuc)
        sub_nuc = [letter for letter in sub_nuc if letter in 'ACTG' and letter != ""]
        if changes[0] == ">" and changes[1] == ">" and changes[2] == "=":
            sub_one, sub_two = int(sub[0]), int(sub[4])
        if changes[0] == ">" and changes[1] == "=" and changes[2] == ">":
            sub_one, sub_two = int(sub[0]), int(sub[6])
            adjacent = False
        if changes[0] == "=" and changes[1] == ">" and changes[2] == ">":
            sub_one, sub_two = int(sub[2]), int(sub[6])
        # get substituted nucleotides
        sub_one_nuc, sub_two_nuc = sub_nuc[1], sub_nuc[3]
        if adjacent:
            standard_hgvs = "c." + str(sub_one) + "_" + str(sub_two) + "delins" + sub_one_nuc + sub_two_nuc
        else:
            standard_hgvs = "c.[" + str(sub_one) + target_seq[sub_one - 1] + ">" + sub_one_nuc + ";" \
                            + str(sub_two) + target_seq[sub_two - 1] + ">" + sub_two_nuc + "]"

    elif count == 3:
        # variant_codon has three nucleotide substitutions
        # get indices of nucleotide substitutions
        sub = re.split("\[", hgvs_legacy)[1]
        sub = re.split("[a-zA-Z.>;=]", sub)
        # get nucleotides of substitutions
        sub_nuc = hgvs_legacy.split("]")[0]
        sub_nuc = re.split("\[", sub_nuc)[1]
        sub_nuc = re.split("[a-z0-9.>;=]", sub_nuc)
        sub_nuc = [letter for letter in sub_nuc if letter in 'ACTG' and letter != ""]
        sub_one, sub_two, sub_three = int(sub[0]), int(sub[4]), int(sub[8])
        # get substituted nucleotides
        sub_one_nuc, sub_two_nuc, sub_three_nuc = sub_nuc[1], sub_nuc[3], sub_nuc[5]
        standard_hgvs = "c." + str(sub_one) + "_" + str(sub_three) + "delins" + sub_one_nuc + sub_two_nuc + sub_three_nuc

    else:
        return None

    return standard_hgvs
