from collections import Counter
import re


def legacy_to_mave_hgvs(hgvs):
    """
    This helper function takes in an hgvs formatted string in c.[1C>A;2=;3=] format and returns the
    indices in the codon that the substitutions occurred as well as the variant nucleotide

    Parameters
    ----------
    hgvs (string): hgvs formatted string

    Returns
    -------
    sub_one, sub_two, sub_three (int): index of nucleotide substitution in codon, None if no substitution
    sub_one_nuc, sub_two_nuc, sub_three_nuc (string): variant nucleotide, None if no substitution
    """

    # determine which bases had a change
    changes = [letter for letter in hgvs if letter in "=>"]
    # count instances of changes
    count = Counter(changes)
    count = count[">"]

    if count == 0:
        # variant_codon is wild-type
        sub_one = None  # no nucleotide substitutions
        sub_two = None
        sub_three = None
        sub_one_nuc = None
        sub_two_nuc = None
        sub_three_nuc = None

    elif count == 1:
        # variant_codon has one nucleotide substitution
        # get indices of nucleotide substitutions
        sub = re.split("\[", hgvs)[1]
        sub = re.split("[a-zA-Z.>;=]", sub)
        # get nucleotides of substitutions
        sub_nuc = hgvs.split("]")[0]
        sub_nuc = re.split("\[", sub_nuc)[1]
        sub_nuc = re.split("[a-z0-9.>;=]", sub_nuc)
        sub_nuc = [letter for letter in sub_nuc if letter in 'ACTG' and letter != ""]
        if changes[0] == ">" and changes[1] == "=" and changes[2] == "=":
            sub_one = int(sub[0])
        if changes[0] == "=" and changes[1] == ">" and changes[2] == "=":
            sub_one = int(sub[2])
        if changes[0] == "=" and changes[1] == "=" and changes[2] == ">":
            sub_one = int(sub[4])
        sub_one = (sub_one % 3) - 1  # index of nucleotide substitution
        sub_one_nuc = sub_nuc[1]
        # set other possible indices for codon substitution to None
        sub_two = None
        sub_two_nuc = None
        sub_three = None
        sub_three_nuc = None

    elif count == 2:
        # variant_codon has two nucleotide substitutions
        # get indices of nucleotide substitutions
        sub = re.split("\[", hgvs)[1]
        sub = re.split("[a-zA-Z.>;=]", sub)
        # get nucleotides of substitutions
        sub_nuc = hgvs.split("]")[0]
        sub_nuc = re.split("\[", sub_nuc)[1]
        sub_nuc = re.split("[a-z0-9.>;=]", sub_nuc)
        sub_nuc = [letter for letter in sub_nuc if letter in 'ACTG' and letter != ""]
        if changes[0] == ">" and changes[1] == ">" and changes[2] == "=":
            sub_one, sub_two = int(sub[0]), int(sub[4])
        if changes[0] == ">" and changes[1] == "=" and changes[2] == ">":
            sub_one, sub_two = int(sub[0]), int(sub[6])
        if changes[0] == "=" and changes[1] == ">" and changes[2] == ">":
            sub_one, sub_two = int(sub[2]), int(sub[6])
        sub_one = (sub_one % 3) - 1  # index of nucleotide substitution
        sub_two = (sub_two % 3) - 1  # index of nucleotide substitution
        sub_one_nuc, sub_two_nuc = sub_nuc[1], sub_nuc[3]
        # set other possible indices for codon substitution to None
        sub_three = None
        sub_three_nuc = None

    elif count == 3:
        # variant_codon has three nucleotide substitutions
        # get indices of nucleotide substitutions
        sub = re.split("\[", hgvs)[1]
        sub = re.split("[a-zA-Z.>;=]", sub)
        # get nucleotides of substitutions
        sub_nuc = hgvs.split("]")[0]
        sub_nuc = re.split("\[", sub_nuc)[1]
        sub_nuc = re.split("[a-z0-9.>;=]", sub_nuc)
        sub_nuc = [letter for letter in sub_nuc if letter in 'ACTG' and letter != ""]
        sub_one, sub_two, sub_three = int(sub[0]), int(sub[4]), int(sub[8])
        sub_one_nuc, sub_two_nuc, sub_three_nuc = sub_nuc[1], sub_nuc[3], sub_nuc[5]
        sub_one = (sub_one % 3) - 1  # index of nucleotide substitution
        sub_two = (sub_two % 3) - 1  # index of nucleotide substitution
        sub_three = (sub_three % 3) - 1  # index of nucleotide substitution

    else:
        return None

    return sub_one, sub_two, sub_three, sub_one_nuc, sub_two_nuc, sub_three_nuc