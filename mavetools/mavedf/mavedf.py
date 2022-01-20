from mavetools.mavedf.df_to_pandas import df_to_pandas

from mavetools.mavedf.legacy_to_mave_new import legacy_to_mave_hgvs_nt
from mavehgvs.variant import Variant
from mavetools.mavedf.mutation_type import *


class MaveDf:
    """
    The MaveDf object consists of a meta_data attribute, a pandas_df attribute, and various functions
    designed to manipulate the pandas_df attribute.
    """

    def __init__(self, df):
        """
        Constructor
        This method instantiates the MaveDf object and assigns values to the
        attributes pandas_df and meta_dict

        Parameters
        ----------
        df : string
            the filename.csv to be converted to pandas df
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
        This method also updates the hgvs from legacy to mave hgvs, if it is not already done.

        Parameters
        __________
        target_seq : string
            target sequence

        Raises
        ______
        TypeError
            if target_seq is not string
        ValueError
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

            # check for legacy hgvs format (i.e., c.[1C>A;2=;3=]) and update if needed
            self.pandas_df["hgvs_nt"][i] = legacy_to_mave_hgvs_nt(
                self.pandas_df["hgvs_nt"][i]
            )

            # new implementation
            # instantiate mavevariant
            # variant = MaveVariant(self.pandas_df["hgvs_nt"][i], target_seq)
            # self.pandas_df["hgvs_nt"][i] = variant.mave_hgvs
            print("variant = " + self.pandas_df["hgvs_nt"][i])

            # get hgvs_nt
            hgvs = self.pandas_df["hgvs_nt"][i]

            # identify variant_position and get codon_number associated with it

            if is_wild_type(hgvs):  # variant_codon is wild-type
                codon_number = None
                target_codon = None
            else:  # any other variant change
                # instantiate Variant object
                variant = Variant(hgvs)
                # get variant position and convert to int
                if type(variant.positions) == list:  # multiple positions values exist
                    variant_position = int(str(variant.positions[0]))
                elif type(variant.positions) == tuple:
                    variant_position = int(str(variant.positions[0]))
                else:  # only one value for positions
                    variant_position = int(str(variant.positions))
                # now that we have the variant_position, get codon_number
                codon_number = round((variant_position / 3) + 0.5)
                # use codon_number to get target_codon from target_seq
                target_codon = target_seq[(codon_number - 1) * 3 : codon_number * 3]

            # determine sequence of variant_codon

            if is_wild_type(hgvs):  # variant_codon is wild-type
                variant_codon = target_codon
                sub_one = None  # no nucleotide substitutions
            elif is_deletion(hgvs):  # target_codon was deleted
                variant_codon = None
                sub_one = None  # no nucleotide substitutions
            elif is_substitution_one_base(
                hgvs
            ):  # variant_codon has one nucleotide substitution
                # instantiate Variant object
                variant = Variant(hgvs)
                # get index of nucleotide substitution
                sub_one = int(str(variant.positions)) % 3 - 1
                # get nucleotide of substitution
                sub_one_nuc = variant.sequence[1]
                # set other possible indices for codon substitution to None
                sub_two = None
                sub_three = None
            elif is_substitution_two_bases_nonadjacent(
                hgvs
            ):  # variant has two nucleotide substitutions, non-adjacent
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
                if (
                    len(sub_nucs) == 2
                ):  # variant codon has two adjacent nucleotide substitutions
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

            # using data generated above (substituted nucleotides and indices in codon), construct variant_codon

            # only assign variant_codon if nucleotide substitution occurred
            if sub_one is not None:
                # declare and initialize variant_codon
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
