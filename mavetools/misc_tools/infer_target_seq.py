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

    return target_seq