import pandas as pd
import numpy as np
import re
from mavetools.new_tools.df_to_pandas import df_to_pandas
from collections import Counter


def add_variant_data(df, target_seq, drop_accession=False, ret_meta=False):
    """
    This function takes in a score or count dataframe and a target sequence and converts coding
    variants into codon changes.

    Parameters
    __________
    df (string): filename.csv file that contains the score or count dataframe
    target_seq (string): target sequence
    drop_accession (bool): True if user wants to drop the accession numbers
        default value = False
    ret_meta (bool): True if user wants to return metadata as dictionary from df
        default value = False

    Returns
    -------
    df_with_variant_data (pandas.df): df with three additional columns:
        target_codon, codon_number, variant_codon
    meta_dict (dictionary): **optional** if ret_meta = True
        dictionary containing formatted metadata (comments) at the top of df

    Raises
    ______
    ValueError: if filename is not in "filename.csv" format
    ValueError: if target_seq is not string of solely characters ACTG
    ValueError: if drop_accession or ret_meta are not bool
    """

    # check for ValueError
    # if df is not in filename.csv format
    if not df.endswith(".csv"):
        raise ValueError("df must be csv file")
    # if target_seq is not string containing solely characters ACTG
    if not isinstance(target_seq, str):
        raise ValueError("target_seq must be string")
    check_chars = [letter in "ACTG" for letter in target_seq]
    if False in check_chars:
        raise ValueError("target_seq is invalid")
    # if drop_accession is not bool
    if not isinstance(drop_accession, bool):
        raise ValueError("drop_accession must be boolean value")
    # if ret_meta is not bool
    if not isinstance(ret_meta, bool):
        raise ValueError("ret_meta must be boolean value")

    # format df in pandas according to arguments
    df_with_variant_data, meta_dict = df_to_pandas(df, drop_accession)

    # add three columns to df
    df_with_variant_data = df_with_variant_data.assign(target_codon="NA")
    df_with_variant_data = df_with_variant_data.assign(codon_number="NA")
    df_with_variant_data = df_with_variant_data.assign(variant_codon="NA")

    # iterate through df and parse hgvs_nt column to get data
    for i in range(len(df_with_variant_data["hgvs_nt"])):

        # identify location and get codon number associated with loc
        codon_location = df_with_variant_data["hgvs_nt"][i]

        # check for alternative hgvs format _______
        # changes will be 3 if in  alternative hgvs format
        changes = [letter for letter in codon_location if letter in "=>"]

        if codon_location.startswith("_wt"):
            # variant_codon is wild-type
            codon_number = None
            target_codon = None
        elif codon_location[2] == "[":
            # if first and third base of codon changed
            # isolate codon location
            codon_location = re.split("\[", codon_location)[1]
            codon_location = int(re.split("[a-zA-Z.>;_=]", codon_location)[0])
            # now that we have codon_location, get codon_number
            codon_number = round((codon_location / 3) + 0.5)
            # use codon_number to get target_codon from target_seq
            target_codon = target_seq[(codon_number - 1) * 3:codon_number * 3]
        else:  # any other variant change
            # isolate codon location
            codon_location = int(re.split("[a-zA-Z.>;_=]", codon_location)[2])
            # now that we have codon_location, get codon_number
            codon_number = round((codon_location/3)+0.5)
            # use codon_number to get target_codon from target_seq
            target_codon = target_seq[(codon_number - 1) * 3:codon_number * 3]

        # determine sequence of variant_codon

        # get hgvs_nt
        hgvs = df_with_variant_data["hgvs_nt"][i]

        if hgvs.startswith("_wt"):
            # variant_codon is wild-type
            variant_codon = target_codon
            sub_one = None  # no nucleotide substitutions

        elif len(changes) == 3:
            # get substitutions for alternative hgvs format
            sub_one, sub_two, sub_three, sub_one_nuc, sub_two_nuc, sub_three_nuc = parse_additional_hgvs_format(hgvs)
            # check for wild-type in alternative hgvs format
            if sub_one is None and sub_two is None and sub_three is None:
                variant_codon = target_codon

        elif hgvs.endswith("del"):
            # target_codon was deleted
            variant_codon = None
            sub_one = None  # no nucleotide substitutions

        elif hgvs[-2] == ">":
            # variant_codon has one nucleotide substitution
            # get index of nucleotide substitution
            sub_one = int(re.split("[a-zA-Z.>;=]", hgvs)[2])  # location of nucleotide in target_seq
            sub_one = (sub_one % 3) - 1 # index of nucleotide substitution
            # get nucleotide of substitution
            sub_one_nuc = hgvs[-1]
            # set other possible indices for codon substitution to None
            sub_two = None
            sub_three = None

        elif hgvs[-1] == "]":  # and len(changes) != 3:  # do not do for alternative hgvs format
            # variant_codon has two nucleotide substitutions, non-adjacent
            # get indices of nucleotide substitutions
            sub = re.split("\[", hgvs)[1]
            sub = re.split("[a-zA-Z.>;=]", sub)
            sub_one, sub_two = int(sub[0]), int(sub[4])  # location of nucleotides in target_seq
            sub_one, sub_two = (sub_one % 3) - 1, (sub_two % 3) - 1  # indices of nucleotide substitutions
            # get nucleotides of substitutions
            sub_nuc = hgvs.split("]")[0]
            sub_nuc = re.split("[a-z0-9.>;=]", sub_nuc)
            sub_one_nuc, sub_two_nuc = sub_nuc[4], sub_nuc[7]
            # set other possible indices for codon substitution to None
            sub_three = None

        else:  # len(changes) != 3:  # any change leftover that is not in alternative format
            # variant_codon has two or three adjacent nucleotide substitutions
            # get index of first codon substitution
            sub_one = int(re.split("[a-zA-Z.>;_=]", hgvs)[2])  # location of first substitution in target_seq
            sub_one = (sub_one % 3) - 1  # index of first nucleotide substitution
            # get string of substituted nucleotides
            sub_nucs = re.split("[a-z0-9.>;_=]", hgvs)[-1]
            if len(sub_nucs) == 2:  # variant codon has two adjacent nucleotide substitutions
                # assign additional nucleotide substitution indices
                sub_two = sub_one + 1
                # get nucleotides of substitutions
                sub_one_nuc = sub_nucs[0]
                sub_two_nuc = sub_nucs[1]
                # set other possible indices for codon substitution to None
                sub_three = None
            else:  # variant has three adjacent nucleotide substitutions
                # assign additional nucleotide substitution indices
                sub_two = sub_one + 1
                sub_three = sub_two + 1
                # get nucleotides of substitutions
                sub_one_nuc = sub_nucs[0]
                sub_two_nuc = sub_nucs[1]
                sub_three_nuc = sub_nucs[2]

        # now that we have the type of change, and stored data for change, get variant_codon
        # but only assign variant_codon if nucleotide substitution occurred
        if sub_one is not None:
            variant_codon = ""

            # set first nucleotide of variant_codon
            if sub_one == 0:
                variant_codon = variant_codon + sub_one_nuc
            else:
                variant_codon = variant_codon + target_codon[0]
            # set second nucleotide of variant_codon
            if sub_one == 1:
                variant_codon = variant_codon + sub_one_nuc
            elif sub_two == 1:
                variant_codon = variant_codon + sub_two_nuc
            else:
                variant_codon = variant_codon + target_codon[1]
            # set third nucleotide of variant_codon
            if sub_one == -1 or sub_one == 2:
                variant_codon = variant_codon + sub_one_nuc
            elif sub_two == -1 or sub_two == 2:
                variant_codon = variant_codon + sub_two_nuc
            elif sub_three == -1 or sub_three == 2:
                variant_codon = variant_codon + sub_three_nuc
            else:
                variant_codon = variant_codon + target_codon[2]

        # add values for target_codon, codon_number, and variant_codon to this row
        df_with_variant_data.at[i, "target_codon"] = target_codon
        df_with_variant_data.at[i, "codon_number"] = codon_number
        df_with_variant_data.at[i, "variant_codon"] = variant_codon

    # return new df_with_variant_data and meta_dict if desired
    if ret_meta:
        return df_with_variant_data, meta_dict
    else:
        return df_with_variant_data


def parse_additional_hgvs_format(hgvs):
    """
    This helper function takes in an hgvs formatted string in _______ format and returns the indeces in
    the codon that the substitutions occurred and the variant nucleotide

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
