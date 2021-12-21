
def delins_to_subs(target, delins, offset=0):
    """
    The code that converts codon changes to mavehgvs variants prefers to define single events,
    e.g., a deletion-insertion of two bases rather than two substitutions. Many users may prefer
    to look at these data as multiple substitutions instead.

    This function takes delins mavehgvs variants that have matching deletion and insertion length
    and outputs the corresponding multi-variant.

    For many variants this will require the user to provide the appropriate target sequence.
    To improve usability, the user is able to provide a longer target sequence and an offset.

    The function should throw an informative error if the delins is not matched or if the target
    sequence doesn't match (e.g. outside of length bounds, or would result in a target-identical
    change that suggests the wrong target sequence was provided).

    Parameters
    ----------
    target : str
        The target sequence that delins is based on.
    delins : str
        A mavehgvs representation of a single deletion-insertion event.
    offset : int
        The offset of index into target sequence. E.g., you have variants from MaveDB that are part of an
        internal domain of a gene. MaveDB's base 1 is actually base 121. In this case, the user can pass
        as target the full-length FASTA and as offset 120.

    Returns
    -------

    Raises
    ______


    """
    return None