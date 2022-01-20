from mavehgvs.variant import Variant


def infer_target_seq(variant_list):
    """
    This function works for both protein and nucleotide data. Input is a list of mavehgvs Variant objects.

    The function handles target identical and indel variants correctly, likely ignoring them in
    most cases (a possible exception being protein data with two adjacent residues).

    Parameters
    ----------
    variant_list : list
        Consists of mavehgvs Variant objects.

    Returns
    -------
    target_seq : str
        The inferred target sequence based on the variant list provided.

    Raises
    ______
    Value Error
        The function outputs an appropriate error or warning message if the variant data has a gap (e.g. if
        there are variants for positions 1-10 and 12-15 but not 11). The function would also throw an error if
        there are conflicting target residues for the same position.
    """
    aa_dict = {"A": "Ala", "C": "Cys", "D": "Asp", "E": "Glu", "F": "Phe",
               "G": "Gly", "H": "His", "I": "Ile", "K": "Lys", "L": "Leu",
               "M": "Met", "N": "Asn", "P": "Pro", "Q": "Gln", "R": "Arg",
               "S": "Ser", "T": "Thr", "V": "Val", "W": "Trp", "Y": "Tyr",
               "*": "Ter", "WTSYN": ""}

    # declare target_seq
    target_seq = ""

    # iterate through variant list
    for variant in variant_list:
        # check if Variaint is not a multivariant and is a substitution
        if variant.is_multi_variant() is False:
            if variant.variant_type == 'sub':
                # check if variant is protein data
                if str(variant.positions)[0] in "ACGHILMOPSTV":  # protein variant
                    # get position, target_base of variant
                    position = int(str(variant.positions)[3:])
                    target_aa = str(variant.positions)[:-len(str(position))]
                    # check if position is beyond current length of target_seq
                    while int(position) > len(target_seq):
                        # append target_seq with N until we reach position
                        target_seq = target_seq + 'N'
                    # now that they are the same length, add target base
                    target_aa = list(aa_dict.keys())[list(aa_dict.values()).index(target_aa)]
                    target_seq = target_seq[0:position - 1] + target_aa + target_seq[position:]
                else:  # nucleotide variant
                    # get position, target_base of variant
                    position = int(str(variant.positions))
                    target_base = str(variant.sequence[0])
                    # check if position is beyond current length of target_seq
                    while int(position) > len(target_seq):
                        # append target_seq with N until we reach position
                        target_seq = target_seq + 'N'
                    # now that they are the same length, add target base
                    target_seq = target_seq[0:position-1] + target_base + target_seq[position:]
            if variant.variant_type == 'delins':
                # do this
                continue
        if variant.is_multi_variant():
            for i in range(len(variant.sequence)):
                if type(variant.sequence[i]) is tuple:
                    if str(variant.positions[i])[0] in "ACGHILMOPSTV":  # protein variant
                        # get position, target_base of variant
                        position = int(str(variant.positions[i])[3:])
                        target_aa = str(variant.positions[i])[:-len(str(position))]
                        # check if position is beyond current length of target_seq
                        while int(position) > len(target_seq):
                            # append target_seq with N until we reach position
                            target_seq = target_seq + 'N'
                        # now that they are the same length, add target base
                        target_aa = list(aa_dict.keys())[list(aa_dict.values()).index(target_aa)]
                        target_seq = target_seq[0:position - 1] + target_aa + target_seq[position:]
                    else:  # nucleotide variant
                        position = int(str(variant.positions[i]))
                        target_base = str(variant.sequence[i][0])
                        # check if position is beyond current length of target_seq
                        while int(position) > len(target_seq):
                            # append target_seq with N until we reach position
                            target_seq = target_seq + 'N'
                        # now that they are the same length, add target base
                        target_seq = target_seq[0:position - 1] + target_base + target_seq[position:]

    if "N" in target_seq:
        raise ValueError("Gaps in target seq.")

    return target_seq
