import re
from mavetools.mavedf.df_to_pandas import df_to_pandas
from mavetools.legacy_to_mave_hgvs.legacy_to_mave_hgvs import legacy_to_mave_hgvs
from mavehgvs.variant import Variant


class MaveDf:

    def __init__(self, df):
        """
        Constructor
        This function instantiates the MaveDf object and assigns values to the
        attributes pandas_df and meta_dict

        Parameters
        ----------
        df (string): the filename.csv to be converted to pandas df
        """
        # convert df to pandas
        self.pandas_df, self.meta_dict = df_to_pandas(df, ret_meta=True)

    def drop_accession(self):
        """
        This method drops the accession numbers from the MaveDf object.
        """
        del self.pandas_df["accession"]

    def add_variant_data(self, target_seq):
        """
        This method takes in a target sequence and converts coding variants into codon changes.
        These changes are stored in three additional MaveDf columns: target_codon, codon_number, and variant_codon.

        Parameters
        __________
        target_seq (string): target sequence

        Raises
        ______
        TypeError:
            if target_seq is not string
        ValueError:
            if target_seq is not made solely of characters ACTG
        """
        # check for TypeError
        # if target_seq is not string
        if not isinstance(target_seq, str):
            raise TypeError("target_seq must be string")

        # check for ValueError
        # if target_seq is not made solely of characters ACTG
        check_chars = [letter in "ACTG" for letter in target_seq]
        if False in check_chars:
            raise ValueError("target_seq is invalid")

        # add three columns to df
        self.pandas_df = self.pandas_df.assign(target_codon="NA")
        self.pandas_df = self.pandas_df.assign(codon_number="NA")
        self.pandas_df = self.pandas_df.assign(variant_codon="NA")

        # iterate through df and parse hgvs_nt column to get data
        for i in range(len(self.pandas_df["hgvs_nt"])):

            # identify location and get codon number associated with loc
            codon_location = self.pandas_df["hgvs_nt"][i]

            # check for alternative hgvs format: c.[1C>A;2=;3=]
            # changes will be 3 if in alternative hgvs format
            changes = [letter for letter in codon_location if letter in "=>"]

            if codon_location.startswith("_wt"):  # variant_codon is wild-type
                codon_number = None
                target_codon = None
            elif codon_location[2] == "[":  # if first and third base of codon changed
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
            hgvs = self.pandas_df["hgvs_nt"][i]

            if hgvs.startswith("_wt"):  # variant_codon is wild-type
                variant_codon = target_codon
                sub_one = None  # no nucleotide substitutions

            elif len(changes) == 3:  # hgvs is in legacy format
                # get substitutions for alternative hgvs format: c.[1C>A;2=;3=]
                sub_one, sub_two, sub_three, sub_one_nuc, sub_two_nuc, sub_three_nuc = legacy_to_mave_hgvs(hgvs)
                # check for wild-type in alternative hgvs format
                if sub_one is None and sub_two is None and sub_three is None:
                    variant_codon = target_codon

            elif hgvs.endswith("del"):  # target_codon was deleted
                variant_codon = None
                sub_one = None  # no nucleotide substitutions

            elif hgvs[-2] == ">":  # variant_codon has one nucleotide substitution
                # instantiate Variant object
                variant = Variant(hgvs)
                # get index of nucleotide substitution
                sub_one = int(str(variant.positions)) % 3 - 1
                # get nucleotide of substitution
                sub_one_nuc = variant.sequence[1]
                # set other possible indices for codon substitution to None
                sub_two = None
                sub_three = None

            elif hgvs[-1] == "]":  # variant_codon has two nucleotide substitutions, non-adjacent
                # instantiate Variant object
                variant = Variant(hgvs)
                # get indices of nucleotide substitutions
                sub_one = int(str(variant.positions[0])) % 3 - 1
                sub_two = int(str(variant.positions[1])) % 3 - 1
                # get nucleotides of substitutions
                sub_one_nuc = variant.sequence[0][1]
                sub_two_nuc = variant.sequence[1][1]
                # set other possible indices for codon substitution to None
                sub_three = None

            else:  # variant_codon has two or three adjacent nucleotide substitutions
                # instantiate Variant object
                variant = Variant(hgvs)
                variant_codon = variant.sequence
                # get index of first codon substitution
                sub_one = int(str(variant.positions[0])) % 3 - 1
                # get string of substituted nucleotides
                sub_nucs = variant.sequence
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
            self.pandas_df.at[i, "target_codon"] = target_codon
            self.pandas_df.at[i, "codon_number"] = codon_number
            self.pandas_df.at[i, "variant_codon"] = variant_codon
