from mavehgvs.variant import Variant

def generate_mavehgvs_string(prefix, wt, mutant, index):
    """
    This function generates a mavehgvs formatted string from a valid prefix, the wt amino acid, codon or base, the
    mutant amino acid, codon or base, and the index. An error will be thrown if incompatible values are passed
    together.

    Parameters
    ----------
    prefix: List[str]
        A list of mavehgvs prefixes describing the variant return types (e.g., "c", "p", et cetera).
    wt: str
        Wildtype amino acid or codon or base.
    mutant: str
        Mutant amino acid, codon or base.
    index: int
        Location of variant in target sequence.

    Returns
    -------
    List[str]
        A list of mavehgvs formatted strings.

    Raises
    ______
    TypeError
        If index is not an int.
    ValueError
        If prefix is invalid.
    ValueError
        If wt or mutant are invalid.
    ValueError
        If incompatible parameters are passed together.
    """
    # validate arguments

    # if p is in prefix


    # if c or n is in prefix

