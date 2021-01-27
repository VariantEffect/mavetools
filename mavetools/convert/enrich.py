from typing import Optional

from mavehgvs import Variant, MaveHgvsParseError
from fqfa import infer_sequence_type, translate_dna
from fqfa.constants import AA_CODES

__all__ = ["seqid_to_variant"]


def seqid_to_variant(seqid: str, wtseq: str) -> Variant:
    """Convert an Enrich seqid string to a mavehgvs Variant object.

    Enrich seqids consist of two comma-delimited lists separated by a '-' character.
    The first list contains 0-indexed positions and the second list contains single-letter amino acid codes.

    By default, Enrich reports amino acid variants.
    This function will accept nucleotide variants as well.

    Parameters
    ----------
    seqid : str
        The Enrich seqid to convert.
    wtseq : str
        The wild-type sequence used in the experiment.
        This can be a nucleotide sequence that will be translated or a protein sequence.

    Returns
    -------
    Variant
        A Variant object describing the substitutions in this seqid.

    Raises
    ------
    ValueError
        If the seqid is not in the expected format.
    ValueError
        If the wild-type sequence is not valid for the given seqid or contains partial codons.

    """
    try:
        positions, alt_amino_acids = (x.split(",") for x in seqid.split("-"))
    except ValueError:
        raise ValueError("incorrectly formatted seqid")

    if len(positions) != len(alt_amino_acids):
        raise ValueError("incorrectly formatted seqid")

    try:
        positions = [int(x) + 1 for x in positions]
        alt_amino_acids = [AA_CODES[x] for x in alt_amino_acids]
    except (ValueError, KeyError):
        raise ValueError("incorrectly formatted seqid")

    wtseq_type = infer_sequence_type(wtseq, report_iupac=True)
    if wtseq_type == "dna":
        wtseq, remainder = translate_dna(wtseq)
        if remainder is not None:
            raise ValueError("wild-type sequence contains incomplete codons")
    elif wtseq_type == "protein":
        pass
    else:
        raise ValueError("invalid wild-type sequence characters")

    if any(x > len(wtseq) for x in positions):
        raise ValueError("wild-type sequence is too short")

    variant_dicts = list()
    for pos, aa in zip(positions, alt_amino_acids):
        variant_dicts.append(
            {
                "variant_type": "sub",
                "prefix": "p",
                "position": pos,
                "target": AA_CODES[wtseq[pos - 1]],
                "variant": aa,
            }
        )

    if len(variant_dicts) == 1:
        return Variant(variant_dicts[0], targetseq=wtseq)
    else:
        return Variant(variant_dicts, targetseq=wtseq)
