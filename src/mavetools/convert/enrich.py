from fqfa.constants import AA_CODES
from mavehgvs import Variant

__all__ = ["seqid_to_variant"]


def seqid_to_variant(seqid: str, wtseq: str) -> Variant:  # noqa: max-complexity 11
    """Convert an Enrich seqid string to a mavehgvs Variant object.

    Enrich seqids consist of two comma-delimited lists separated by a '-' character.
    The first list contains 0-indexed positions and the second list contains single-letter amino acid codes.

    Parameters
    ----------
    seqid : str
        The Enrich seqid to convert.
    wtseq : str
        The wild-type sequence used in the experiment.
        This must be a protein sequence.

    Returns
    -------
    Variant
        A Variant object describing the substitutions in this seqid.

    Raises
    ------
    ValueError
        If the seqid is not in the expected format.
    ValueError
        If the wild-type sequence is not valid for the given seqid.

    """
    try:
        positions, alt_amino_acids = (x.split(",") for x in seqid.split("-"))
    except ValueError:
        raise ValueError("unable to parse seqid")

    if len(positions) != len(alt_amino_acids):
        raise ValueError("mismatched number of positions and amino acids in seqid")

    try:
        positions = [int(x) + 1 for x in positions]
        alt_amino_acids = [AA_CODES[x] for x in alt_amino_acids]
    except ValueError:
        raise ValueError("unable to parse positions in seqid")
    except KeyError:
        raise ValueError("unexpected amino acid symbol in seqid")

    if max(positions) > len(wtseq):
        raise ValueError("wild-type sequence is too short")

    variant_dicts = list()
    for pos, aa in zip(positions, alt_amino_acids):
        wt_aa = AA_CODES[wtseq[pos - 1]]
        if wt_aa == aa:
            raise ValueError("synonymous changes not expected in seqid")
        variant_dicts.append(
            {
                "variant_type": "sub",
                "prefix": "p",
                "position": pos,
                "target": wt_aa,
                "variant": aa,
            }
        )

    if len(variant_dicts) == 1:
        return Variant(variant_dicts[0], targetseq=wtseq)
    else:
        return Variant(variant_dicts, targetseq=wtseq)
