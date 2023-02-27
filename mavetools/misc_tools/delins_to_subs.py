
def delins_to_subs(target, delins, offset=0):
    """
    The code that converts codon changes to mavehgvs variants prefers to define single events
    (e.g., a deletion-insertion of two bases rather than two substitutions). Many users may prefer
    to look at these data as multiple substitutions instead.

    This function takes delins mavehgvs variants that have matching deletion and insertion length
    and outputs the corresponding multi-variant.

    For many variants this will require the user to provide the appropriate target sequence.
    To improve usability, the user is able to provide a longer target sequence and an offset.

    Parameters
    ----------
    target : str
        The target sequence that delins is based on.
    delins : str
        A mavehgvs representation of a single deletion-insertion event.
    offset : int
        The offset of index into target sequence (e.g., you have variants from MaveDB that are part of an
        internal domain of a gene. MaveDB's base 1 is actually base 121. In this case, the user can pass
        as target the full-length FASTA and as offset 120).

    Returns
    -------
    multi_variant : str
        The corresponding multi-variant based on the delins single event

    Raises
    ______
    ValueError
        If the target sequence is identical to the variant sequence.
    ValueError
        If the variant sequence is outside the length of bounds of the target sequence.
    """
    # save first part of delins string and append start of return string
    multi_variant = delins[0:2] + "["

    # find range of indices in delins
    start_delins, end_delins = delins.split("_")
    start_delins = int(''.join(i for i in start_delins if i.isdigit()))
    end_delins = int(''.join(i for i in end_delins if i.isdigit()))

    # throw error if the end of variant is outside the bounds of the target
    if end_delins + offset > len(target):
        raise ValueError("variant outside of length of bounds of target")

    # save bases in delins
    bases_delins = ''.join(i for i in delins if i.isupper())

    # with range, get appropriate subsequence from target
    bases_target = target[start_delins-1+offset:end_delins+offset]

    # throw error if bases_delins and bases_target are identical
    if bases_delins == bases_target:
        raise ValueError("variant sequence matches target sequence")

    # compare subsequence in target with bases in delins
    for i in range(len(bases_target)):
        # where they differ, add to return string, with appropriate index
        if bases_target[i] != bases_delins[i]:
            multi_variant = multi_variant + str(start_delins+i) + bases_target[i] + ">" + bases_delins[i] + ";"

    # close off multi_variant
    multi_variant = multi_variant[:-1] + "]"

    # return multi_variant
    return multi_variant