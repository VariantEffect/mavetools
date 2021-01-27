from typing import Optional
from fqfa.validator import dna_bases_validator

__all__ = ["codon_sub_to_mavehgvs"]


def codon_sub_to_mavehgvs(
    target_codon: str,
    variant_codon: str,
    aa_position: int,
    target_id: Optional[str] = None,
) -> str:
    """Create a MAVE-HGVS coding variant string describing the change between two codons.

    Parameters
    ----------
    target_codon : str
        Three-letter string for the target sequence codon.
    variant_codon : str
        Three-letter string for the variant codon.
    aa_position : int
        The amino acid position for this codon in the target.
        This will be used to calculate the nucleotide positions.
    target_id : Optional[str]
        Optional target identifier for the resulting variant.

    Returns
    -------
    str
        MAVE-HGVS string for the substitution described by this codon.

    Raises
    ------
    ValueError
        If either codon is not valid.

    """
    if len(target_codon) != 3 or not dna_bases_validator(target_codon):
        raise ValueError("invalid target codon")
    if len(variant_codon) != 3 or not dna_bases_validator(variant_codon):
        raise ValueError("invalid variant codon")

    if target_codon == variant_codon:
        variant_string = "c.="
    else:
        variant_pos = (
            aa_position - 1
        ) * 3 + 1  # nucleotide position of the start of the codon
        changes = [x != y for x, y in zip(target_codon, variant_codon)]
        if sum(changes) == 1:  # single substitution
            codon_pos = changes.index(True)
            variant_string = f"c.{variant_pos + codon_pos}{target_codon[codon_pos]}>{variant_codon[codon_pos]}"
        elif sum(changes) == 2:
            if changes[0] and changes[2]:  # two separate substitutions
                variant_string = (
                    f"c.[{variant_pos}{target_codon[0]}>{variant_codon[0]};"
                    f"{variant_pos + 2}{target_codon[2]}>{variant_codon[2]}]"
                )
            elif changes[0]:  # delins of first two bases
                variant_string = (
                    f"c.{variant_pos}_{variant_pos + 1}delins{variant_codon[:2]}"
                )
            else:  # delins of last two bases
                variant_string = (
                    f"c.{variant_pos + 1}_{variant_pos + 2}delins{variant_codon[1:]}"
                )
        elif sum(changes) == 3:  # full codon delins
            variant_string = f"c.{variant_pos}_{variant_pos + 2}delins{variant_codon}"
        else:  # pragma: nocover
            raise ValueError("invalid codon substitution")

    if target_id is not None:
        return f"{target_id}:{variant_string}"
    else:
        return variant_string
