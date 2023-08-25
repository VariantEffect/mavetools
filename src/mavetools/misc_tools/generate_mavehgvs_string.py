from mavehgvs.variant import Variant
import mavehgvs
import pandas as pd


# TODO: delete this file probably?
def generate_mavehgvs_string(prefix, wt, mutant, position: int, multi_variant=False):
    """
    This function generates a mavehgvs formatted string from a valid prefix, the wt amino acid, codon or base, the
    mutant amino acid, codon or base, and the position. An error will be thrown if incompatible values are passed
    together.

    Parameters
    ----------
    prefix: str
        A mavehgvs prefix describing the variant return type (e.g., "c", "p", et cetera).
    wt: str
        Wildtype amino acid or codon or base.
    mutant: str
        Mutant amino acid, codon or base.
    position: int
        Location of variant in target sequence.

    Returns
    -------
    str
        A mavehgvs formatted string.

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
    if prefix == "p":
        hgvs = protein_mavehgvs(wt, mutant, position)
    # if c or n is in prefix
    elif prefix == "c" or prefix == "n":
        pass  # TODO: this is broken
    return hgvs


def protein_mavehgvs(wt, mutant, position: int):
    if pd.isna(position):
        return None
    if wt == mutant:
        hgvs = "p." + wt + str(position) + "="
    else:
        hgvs = "p." + wt + str(position) + mutant
    # validate variant
    hgvs_validate = mavehgvs.Variant(hgvs)
    return hgvs  # , hgvs_validate
