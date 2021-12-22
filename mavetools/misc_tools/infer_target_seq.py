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

    # declare target_seq
    target_seq = ""

    # iterate through variant list
    for variant in variant_list:
        # check if Variaint is not a multivariant and is a substitution
        if variant.is_multi_variant() is False and variant.variant_type == 'sub':
            # get position, target_base of variant
            position = int(str(variant.positions))
            target_base = str(variant.sequence[0])
            # check if position is beyond current length of target_seq
            #print(str(position))
            while int(position) > len(target_seq):
                # append target_seq with N until we reach position
                target_seq = target_seq + 'N'
            # now that they are the same length, add target base
            #target_seq[position-1] = target_base
            target_seq = target_seq[0:position-1] + target_base + target_seq[position:]



        # check if it is delins or not

        # if it is delins, we cannot infer target seq
        # if it is not, we will be able to, we just want to know if variant type = sub
            #

    return target_seq