import math
from collections import Counter
import re


class MaveVariant:

    def __init__(self, legacy_hgvs, target_seq):
        """
        Constructor
        This method instantiates the MaveVariant object

        Parameters
        ----------
        legacy_mave (string): the legacy formatted hgvs string to be converted
        target_seq (string): target sequence
        """
        # save legacy_hgvs and target_seq as itself in constructor
        self.legacy_hgvs = legacy_hgvs
        self.target_seq = target_seq

        # make dictionary of lists of tuples
        # key = codon number and
        # value = list of (location_in_codon, old base, new base)
        self.mave_hgvs_dict = dict()

        # default value for mave_hgvs is legacy_hgvs, changes when conversion is made
        self.mave_hgvs = None

        self.legacy_to_mave_hgvs_nt()

    def legacy_to_mave_hgvs_nt(self):
        """
        This function converts a legacy hgvs_nt formatted string (i.e., c.[1C>A;2=;3=]) and converts it
        to the standard format (i.e., c.1delinsA). If string is already in standard format, self.mave_hgvs
        is saved as self.legacy_hgvs.

        Parameters
        ----------
        self: MaveHgvsVariant object

        Returns
        -------
        0 (int) if no conversion was necessary
        1 (int) if legacy_hgvs was converted
        """
        # first check for _wt
        if self.legacy_hgvs == "_wt":
            self.mave_hgvs = self.legacy_hgvs
            print("_wt : " + self.mave_hgvs)
            return 0

        # if it doesnt start with a bracket then it is in the correct format
        if self.legacy_hgvs[2] != "[":
            # unless there is an = sign, could be _wt
            if "=" not in self.legacy_hgvs:
                self.mave_hgvs = "_wt"
                return 1
            else:  # no conversion necessary
                self.mave_hgvs = self.legacy_hgvs
                print("in correct format : " + self.mave_hgvs)
                return 0
        else:  # remove everything before and after bracket, including brackets
            # save and keep leading two characters
            leading_chars = self.legacy_hgvs[0:2]
            # remove characters upto first bracket and remove last bracket
            substitutions = self.legacy_hgvs[3:-1]

        # split at semi colon
        # do check for wild-type seq before and after
        count = Counter(substitutions)
        count = count["="]
        substitutions = re.split(";", substitutions)
        # check for wild-type seq
        if count == 3 and len(substitutions) == 3:
            self.mave_hgvs = "_wt"
            return 1

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
                self.mave_hgvs = self.legacy_hgvs
                return 0

            # get substitution and location in codon
            # get locations of nucleotide substitutions
            sub_loc = int(re.split("[a-zA-Z>]", substitution)[0])
            # get nucleotides of substitutions
            sub_nuc = re.split("[a-z0-9>]", substitution)
            sub_nuc = [letter for letter in sub_nuc if letter in 'ACTG' and letter != ""]

            # now that we have the sub_loc, get codon_number
            codon_number = math.ceil(sub_loc / 3)

            # add to mave_hgvs_dict by codon number
            # key = codon number and
            # value = list of [location_in_codon, old_base, new_base]
            key = int(codon_number)
            #location_in_codon = int(sub_loc) % 3 - 1
            location_in_codon = sub_loc
            old_base, new_base = sub_nuc[0], sub_nuc[1]
            value = (location_in_codon, old_base, new_base)
            value_list = [(location_in_codon, old_base, new_base)]
            # update dictionary
            if key not in self.mave_hgvs_dict:  # then add key as empty list
                self.mave_hgvs_dict.update({key: []})
            # add values to list
            # get current list
            current_list = list(self.mave_hgvs_dict.get(key))
            # append list
            current_list.append(value)
            # update dictionary
            self.mave_hgvs_dict.update({key: current_list})

            #if key in self.mave_hgvs_dict:  # first check if key exists
                # if it does, update the value
            #    value_list.insert(0, self.mave_hgvs_dict.get(key))
            # add correct key : value pair to dictionary
            #self.mave_hgvs_dict.update({key: value_list})
            #print(self.mave_hgvs_dict)
        # start constructing string from leading chars saved earlier
        constructing_mave_hgvs = leading_chars + "["
        #print(self.mave_hgvs_dict)
        # loop through dictionary in sorted order and continue constructing mave_hgvs string
        for key in sorted(self.mave_hgvs_dict):  # look at all values for each codon number
            # get the value (list of tuples) for this key (codon number)
            previous_value = None
            current_value = None
            adjacent_values = []  # [(start, end, value)]
            in_construction = False  # keep track of when constructing delins
            for value in self.mave_hgvs_dict[key]:
                #print("VALUE : " + str(value))
                if len(adjacent_values) == 0:
                    #print(value)
                    # add value to adjacent values
                    if type(value) == tuple:
                        adjacent_values.append(value)
                    else:  # value is in list form
                        adjacent_values.append(value[0])
                    #print(adjacent_values)
                else:  # check if next value is adjacent to previous
                    #print(value)
                    if adjacent_values[-1][0] - value[0] == -1:  # value is adjacent to previous value
                        #print(value)
                        # add value to adjacent values
                        adjacent_values.append(value)
                    else:  # value is not adjacent to previous values
                        # add previous values to constructed string
                        if len(adjacent_values) == 1:
                            constructing_mave_hgvs = constructing_mave_hgvs + \
                                                     str(adjacent_values[0][0]) + \
                                                     adjacent_values[0][1] + \
                                                     ">" + \
                                                     adjacent_values[0][2] + \
                                                     ";"
                        else:  # more than one value in adjacent_values
                            constructing_mave_hgvs = constructing_mave_hgvs + \
                                                     str(adjacent_values[0][0]) + \
                                                     "_" + \
                                                     str(adjacent_values[-1][0]) + \
                                                     "delins"
                            for adj_val in adjacent_values:
                                # add new base
                                constructing_mave_hgvs = constructing_mave_hgvs + adj_val[2] + ";"
                        # set previous values as empty
                        adjacent_values.clear()
                        # add new value to adjacent values
                        adjacent_values.append(value)

            # now we have gone through all values, add last values remaining in adjacent values to string
            #print("len " + str(len(adjacent_values)))
            if len(adjacent_values) == 1:
                constructing_mave_hgvs = constructing_mave_hgvs + \
                                         str(adjacent_values[0][0]) + \
                                         adjacent_values[0][1] + \
                                         ">" + \
                                         adjacent_values[0][2] + \
                                         ";"
            else:  # more than one value in adjacent_values
                #print(adjacent_values)
                #print(adjacent_values[0])
                constructing_mave_hgvs = constructing_mave_hgvs + \
                                         str(adjacent_values[0][0]) + \
                                         "_" + \
                                         str(adjacent_values[-1][0]) + \
                                         "delins"
                for adj_val in adjacent_values:
                    # add new base
                    constructing_mave_hgvs = constructing_mave_hgvs + adj_val[2]
                constructing_mave_hgvs = constructing_mave_hgvs + ";"
            # be sure to group adjacent changes in same codon together
            # should I partition by codon number?

        # add closing bracket
        #constructing_mave_hgvs[-1] = "]"
        # string consructed, assign to mave_hgvs
        self.mave_hgvs = constructing_mave_hgvs[:-1] + "]"
        # make sure brakcets are only there when needed
        count = Counter(self.mave_hgvs)
        count = count[";"]
        if count == 0:  # remove brackets
            self.mave_hgvs = leading_chars + self.mave_hgvs[3:-1]
        #print(self.mave_hgvs)
        return 1

    def validate_target_seq(self):
        return None