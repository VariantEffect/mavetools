import math
from collections import Counter
import re


def legacy_to_mave_hgvs_nt(legacy_hgvs):
    """
    This function converts a legacy hgvs_nt formatted string (i.e., c.[1C>A;2=;3=]) and converts it
    to the standard format (i.e., c.1delinsA). If string is already in standard format, the input
    string is returned as-is.

    Parameters
    ----------
    hgvs_legacy : string
        legacy format hgvs string

    Returns
    -------
    mave_hgvs : string
        standard format hgvs string or
        empty string if input string is in invalid format
    """
    # dictionary to keep track of changes
    mave_hgvs_dict = dict()

    # first check for _wt
    if legacy_hgvs == "_wt":
        mave_hgvs = legacy_hgvs
        return mave_hgvs

    # if it doesnt start with a bracket then it is in the correct format
    if legacy_hgvs[2] != "[":
        # unless there is an = sign, could be _wt
        if "=" in legacy_hgvs:
            mave_hgvs = "_wt"
            return mave_hgvs
        else:  # no conversion necessary
            mave_hgvs = legacy_hgvs
            return mave_hgvs
    else:  # remove everything before and after bracket, including brackets
        # save and keep leading two characters
        leading_chars = legacy_hgvs[0:2]
        # remove characters upto first bracket and remove last bracket
        substitutions = legacy_hgvs[3:-1]

    # split at semi colon
    # do check for wild-type seq before and after
    count = Counter(substitutions)
    count = count["="]
    substitutions = re.split(";", substitutions)
    # check for wild-type seq
    if count == 3 and len(substitutions) == 3:
        mave_hgvs = "_wt"
        return mave_hgvs

    # save each string location and base change, if there was a change, i.e not =
    for substitution in substitutions:
        # do the following only if there was a change, not if there is an =
        if "=" in substitution:
            # no conversion necessary
            continue
        # if the substitution is already in mave_hgvs  format, i.e., had del or delins already in it
        # make sure to capture all cases!
        # return 0 and save mave_hgvs and legacy_hgvs
        if "delins" in substitution:
            mave_hgvs = legacy_hgvs
            return mave_hgvs

        # get substitution and location in codon
        # get locations of nucleotide substitutions
        sub_loc = int(re.split("[a-zA-Z>]", substitution)[0])
        # get nucleotides of substitutions
        sub_nuc = re.split("[a-z0-9>]", substitution)
        sub_nuc = [letter for letter in sub_nuc if letter in "ACTG" and letter != ""]

        # now that we have the sub_loc, get codon_number
        codon_number = math.ceil(sub_loc / 3)

        # add to mave_hgvs_dict by codon number
        # key = codon number and
        # value = list of [location_in_codon, old_base, new_base]
        key = int(codon_number)
        # location_in_codon = int(sub_loc) % 3 - 1
        location_in_codon = sub_loc
        old_base, new_base = sub_nuc[0], sub_nuc[1]
        value = (location_in_codon, old_base, new_base)
        value_list = [(location_in_codon, old_base, new_base)]
        # update dictionary
        if key not in mave_hgvs_dict:  # then add key as empty list
            mave_hgvs_dict.update({key: []})
        # add values to list
        # get current list
        current_list = list(mave_hgvs_dict.get(key))
        # append list
        current_list.append(value)
        # update dictionary
        mave_hgvs_dict.update({key: current_list})

    # start constructing string from leading chars saved earlier
    constructing_mave_hgvs = leading_chars + "["
    # loop through dictionary in sorted order and continue constructing mave_hgvs string
    for key in sorted(mave_hgvs_dict):  # look at all values for each codon number
        # get the value (list of tuples) for this key (codon number)
        previous_value = None
        current_value = None
        adjacent_values = []  # [(start, end, value)]
        in_construction = False  # keep track of when constructing delins
        for value in mave_hgvs_dict[key]:
            if len(adjacent_values) == 0:
                # add value to adjacent values
                if type(value) == tuple:
                    adjacent_values.append(value)
                else:  # value is in list form
                    adjacent_values.append(value[0])
            else:  # check if next value is adjacent to previous
                if (
                    adjacent_values[-1][0] - value[0] == -1
                ):  # value is adjacent to previous value
                    # add value to adjacent values
                    adjacent_values.append(value)
                else:  # value is not adjacent to previous values
                    # add previous values to constructed string
                    if len(adjacent_values) == 1:
                        constructing_mave_hgvs = (
                            constructing_mave_hgvs
                            + str(adjacent_values[0][0])
                            + adjacent_values[0][1]
                            + ">"
                            + adjacent_values[0][2]
                            + ";"
                        )
                    else:  # more than one value in adjacent_values
                        constructing_mave_hgvs = (
                            constructing_mave_hgvs
                            + str(adjacent_values[0][0])
                            + "_"
                            + str(adjacent_values[-1][0])
                            + "delins"
                        )
                        for adj_val in adjacent_values:
                            # add new base
                            constructing_mave_hgvs = (
                                constructing_mave_hgvs + adj_val[2] + ";"
                            )
                    # set previous values as empty
                    adjacent_values.clear()
                    # add new value to adjacent values
                    adjacent_values.append(value)

        # now we have gone through all values, add last values remaining in adjacent values to string
        if len(adjacent_values) == 1:
            constructing_mave_hgvs = (
                constructing_mave_hgvs
                + str(adjacent_values[0][0])
                + adjacent_values[0][1]
                + ">"
                + adjacent_values[0][2]
                + ";"
            )
        else:  # more than one value in adjacent_values
            constructing_mave_hgvs = (
                constructing_mave_hgvs
                + str(adjacent_values[0][0])
                + "_"
                + str(adjacent_values[-1][0])
                + "delins"
            )
            for adj_val in adjacent_values:
                # add new base
                constructing_mave_hgvs = constructing_mave_hgvs + adj_val[2]
            constructing_mave_hgvs = constructing_mave_hgvs + ";"
        # be sure to group adjacent changes in same codon together
        # should I partition by codon number?

    # add closing bracket
    # constructing_mave_hgvs[-1] = "]"
    # string consructed, assign to mave_hgvs
    mave_hgvs = constructing_mave_hgvs[:-1] + "]"
    # make sure brackets are only there when needed
    count = Counter(mave_hgvs)
    count = count[";"]
    if count == 0:  # remove brackets
        mave_hgvs = leading_chars + mave_hgvs[3:-1]

    return mave_hgvs
