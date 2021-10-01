def is_wild_type(hgvs):
    """
    This function takes an hgvs formatted string and returns True if the hgvs string indicates
    there was no change from the target sequence.

    Parameters
    ----------
    hgvs : string
        hgvs formatted string

    Returns
    -------
    wt : bool
        True if hgvs string indicates wild type
    """
    wt = False
    if hgvs.startswith("_wt"):
        wt = True
    return wt


def is_deletion(hgvs):
    """
    This function takes an hgvs formatted string and returns True if the hgvs string indicates
    there was a deletion.

    Parameters
    ----------
    hgvs : string
        hgvs formatted string

    Returns
    -------
    deletion : bool
        True if hgvs string is indicates a deletion
    """
    deletion = False
    if hgvs.endswith("del"):
        deletion = True
    return deletion


def is_substitution_one_base(hgvs):
    """
    This function takes an hgvs formatted string and returns True if the hgvs string indicates
    there was a substitution at one base of the codon.

    Parameters
    ----------
    hgvs : string
        hgvs formatted string

    Returns
    -------
    sub_one : bool
        True if hgvs string is indicates a substitution at one base of codon
    """
    sub_one = False
    if hgvs[-2] == ">":
        sub_one = True
    return sub_one


def is_substitution_two_bases_nonadjacent(hgvs):
    """
    This function takes an hgvs formatted string and returns True if the hgvs string indicates
    there were substitutions (non-adjacent) in the codon.

    Parameters
    ----------
    hgvs : string
        hgvs formatted string

    Returns
    -------
    sub_two : bool
        True if hgvs string is indicates a substitution at one base of codon
    """
    sub_two = False
    if hgvs[-1] == "]":
        sub_two = True
    return sub_two
