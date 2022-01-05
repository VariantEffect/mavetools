import re
from collections import OrderedDict

from mavehgvs.patterns import rna, dna, protein, single_variant_re, multi_variant_re

import numpy as np
import pandas as pd
from fqfa.constants.translation.table import CODON_TABLE
from fqfa.constants.iupac.protein import AA_CODES

from . import constants, exceptions


def slicer(seq, size):
    """Slices a string into chunks of `size`."""
    return (seq[pos : pos + size] for pos in range(0, len(seq), size))


def translate_dna(wt_sequence, offset=0):
    """
    Translates a DNA wild-type sequence starting from an `offset`.
    Parameters
    ----------
    wt_sequence : str
        The wild-type sequence to translate.
    offset : int
        Number of bases at the beginning of `wt_sequence` to ignore before
        beginning translation.
    Returns
    -------
    `str`
        The translated wild-type sequence.
    """
    if offset < 0:
        raise ValueError("Offset must not be negative.")

    coding_region = wt_sequence[offset:]
    if len(coding_region) % 3 != 0:
        raise ValueError(
            "Length of sequence derived using an offset of {} is "
            "not a multiple of 3 and cannot be translated.".format(offset)
        )
    protein_seq = ""
    for codon in slicer(coding_region, 3):
        # Let if fail loudly for now
        protein_seq += CODON_TABLE[codon.upper()]
    return protein_seq


def is_null(value):
    """
    Returns `True` if `value` is null, undefined, none, na, n/a, nan or empty.
    """
    value = str(value).strip().lower()
    return (not value) or constants.null_value_re.fullmatch(value) is not None


def format_column(values, astype=float):
    """
    Formats a list of values by replacing null float/int values with
    `np.NaN` or null object values with None. All other values
    are typecast to `astype`.
    Parameters
    ----------
    values : list[Union[float, int]]
        List of values to format.
    astype : callable, optional
        Type-casting callback accepting a single argument.
    Returns
    -------
    list[Any]
        List of values with type returned by `astype` and null values
        replaced with `np.NaN`.
    """
    cast_to_numeric = is_numeric(astype)
    none_type = np.NaN if cast_to_numeric else None
    return [none_type if is_null(v) else astype(v) for v in values]


def is_numeric(dtype):
    """
    Returns `True` if a dtype is a subtype of a `float` or `int`.
    Parameters
    ----------
    dtype: type
        Type to check
    Returns
    -------
    bool
    """
    return np.issubdtype(dtype, np.floating) or np.issubdtype(dtype, np.signedinteger)


class NucleotideSubstitutionEvent(object):
    """
    Parses a nucleotide HGVS_ string into a python class. Can only accept
    basic strings of the format `<prefix>.<position><ref>><alt>`
    Attributes
    ----------
    position : int
        Position of the substitution event.
    ref : str, optional.
        The reference base. May be `None` if the substitution is silent.
    alt : str, optional.
        The mutant base. May be `None` if the substitution is silent.
    silent : bool
        `True` if `ref == alt` or the event has the form `<prefix>.<position>=`
    prefix : str
        Prefix of the variant.
    """

    def __init__(self, variant):
        self.variant = variant.strip()
        match_dna = dna.substitution_re.fullmatch(self.variant)
        match_rna = rna.substitution_re.fullmatch(self.variant)
        if not (match_dna or match_rna):
            raise exceptions.InvalidVariantType(
                "'{}' is not a valid DNA/RNA "
                "substitution event.".format(self.variant)
            )

        match = match_dna if match_dna is not None else match_rna
        self.dict = match.groupdict()
        self.position = int(match.groupdict()[constants.hgvsp_nt_pos])
        self.ref = match.groupdict()[constants.hgvsp_nt_ref]
        self.alt = match.groupdict()[constants.hgvsp_nt_alt]
        self.silent = match.groupdict()[constants.hgvsp_silent] == "="
        self.prefix = variant[0].lower()

        if self.dict.get("utr", None) == "-":
            self.position *= -1

        if match_rna and self.position < 0:
            raise IndexError("RNA positions cannot be negative.")

        if self.ref:
            if match_dna:
                self.ref = self.ref.upper()
        if self.alt:
            if match_dna:
                self.alt = self.alt.upper()
        if self.ref == self.alt:
            self.silent = True

    def __repr__(self):
        return self.format

    @property
    def format(self):
        return "{}.{}".format(self.prefix, self.event)

    @property
    def event(self):
        if self.silent:
            return "{pos}=".format(ref=self.ref, pos=self.position)
        return "{pos}{ref}>{alt}".format(ref=self.ref, pos=self.position, alt=self.alt)

    def codon_position(self, one_based=True):
        """
        Returns the 1-based codon position of a variant.
        Parameters
        ----------
        one_based : bool
            Set as `True` if the variant position expressed in
            1-based coordinates.
        Returns
        -------
        int
        """
        if self.position < 0:
            raise ValueError("Cannot infer codon from negative position.")
        return (self.position - int(one_based)) // 3 + 1

    def codon_frame_position(self, one_based=True):
        """
        Returns the 1-based position of this variant within it's codon.
        Parameters
        ----------
        one_based : bool
            Set as `True` if the variant position expressed in
            1-based coordinates.
        Returns
        -------
        int
        """
        if self.position < 0:
            raise ValueError("Cannot infer codon frame from negative position.")
        return (
            self.position
            - 3 * (self.codon_position(one_based) - 1)
            + int(not one_based)
        )


class ProteinSubstitutionEvent(object):
    """
    Parses a protein HGVS_ string into a python class. Can only accept
    basic strings of the format `p.<ref><position><alt>`
    Attributes
    ----------
    position : int
        Position of the substitution event.
    ref : str, optional.
        The reference amino acid in three-letter-code format.
    alt : str, optional.
        The mutant amino acid in three-letter-code format.
    silent : bool
        `True` if `ref == alt` or the event has the form
        `<prefix>.<ref><position>=`
    prefix : str
        Prefix of the variant.
    """

    def __init__(self, variant):
        self.variant = variant.strip()
        match = protein.substitution_re.fullmatch(self.variant)
        if not match:
            raise exceptions.InvalidVariantType(
                "'{}' is not a valid amino acid "
                "substitution event.".format(self.variant)
            )

        self.dict = match.groupdict()
        self._position = None
        self.position = int(match.groupdict()[constants.hgvsp_pro_pos])
        self.ref = match.groupdict()[constants.hgvsp_pro_ref]
        self.alt = match.groupdict()[constants.hgvsp_pro_alt]
        self.silent = match.groupdict()[constants.hgvsp_silent] == "="
        self.prefix = "p"

        # Normalize to three letter codes
        if self.ref and len(self.ref) == 1:
            self.ref = AA_CODES[self.ref]
        if self.alt and len(self.alt) == 1:
            if self.alt == "?":
                self.alt = "???"
            else:
                self.alt = AA_CODES[self.alt]

        if self.ref and self.silent:
            self.alt = self.ref

    def __repr__(self):
        return self.format

    @property
    def position(self):
        return self._position

    @position.setter
    def position(self, value):
        if value < 1:
            raise ValueError(
                "Protein position cannot be less "
                "Attempted to set {} in {}".format(value, self.variant)
            )
        self._position = value

    @property
    def format(self):
        return "{}.{}".format(self.prefix, self.event)

    @property
    def event(self):
        if self.silent:
            return "{ref}{pos}=".format(ref=self.ref, pos=self.position)
        return "{ref}{pos}{alt}".format(ref=self.ref, pos=self.position, alt=self.alt)


def split_variant(variant):
    """
    Splits a multi-variant `HGVS` string into a list of single variants. If
    a single variant string is provided, it is returned as a singular `list`.
    Parameters
    ----------
    variant : str
        A valid single or multi-variant `HGVS` string.
    Returns
    -------
    list[str]
        A list of single `HGVS` strings.
    """
    prefix = variant[0]
    if len(variant.split(";")) > 1:
        return ["{}.{}".format(prefix, e.strip()) for e in variant[3:-1].split(";")]
    return [variant]


def normalize_variant(variant):
    """
    Replaces `???` for `Xaa` in protein variants and `X` for `N` in
    nucleotide variants to be compliant with the `hgvs` biocommons package.
    Use for enrich and enrich2 inputs.
    Parameters
    ----------
    variant : str, optional.
        HGVS_ formatted string.
    Returns
    -------
    str
    """
    if variant is None:
        return variant
    if (
        protein.single_variant_re.fullmatch(variant)
        or protein.multi_variant_re.fullmatch(variant)
        or protein.predicted_variant_re.fullmatch(variant)
        or protein.any_event_re.fullmatch(variant)
    ):
        # Sub groups of three first.
        variant = re.sub(r"\?{3}", "Xaa", variant)
        # Sub singular next.
        variant = re.sub(r"\?", "X", variant)
    elif (
        dna.single_variant_re.fullmatch(variant)
        or dna.multi_variant_re.fullmatch(variant)
        or dna.any_event_re.fullmatch(variant)
    ):
        variant = variant.replace(r"X", "N")
    elif (
        rna.single_variant_re.fullmatch(variant)
        or rna.multi_variant_re.fullmatch(variant)
        or rna.any_event_re.fullmatch(variant)
    ):
        variant = variant.replace(r"x", "n")
    return variant.strip()


def format_variant(variant):
    """
    Return None for null variant and strips trailing whitespaces.
    Parameters
    ----------
    variant : str, optional.
        HGVS_ formatted string.
    Returns
    -------
    str
    """
    if variant is None:
        return variant
    return variant.strip()


def hgvs_pro_from_event_list(events):
    """
    Convert a list of protein variant events into a single HGVS string. Removes
    duplicates from `events`.
    """
    events = list(OrderedDict.fromkeys([format_variant(e) for e in events]).keys())
    if len(events) == 1:
        mave_hgvs = "p.{}".format(format_variant(events[0]))
    else:
        mave_hgvs = "p.[{}]".format(";".join(events))

    match = protein.single_variant_re.fullmatch(
        mave_hgvs
    ) or protein.multi_variant_re.fullmatch(mave_hgvs)
    if not match:
        raise exceptions.HGVSMatchError(
            "Could not validate parsed variant '{variant}'.".format(variant=mave_hgvs)
        )
    return mave_hgvs


def hgvs_nt_from_event_list(events, prefix):
    """Convert a list of variant events into a single HGVS string."""
    if len(events) == 1:
        mave_hgvs = "{}.{}".format(prefix, format_variant(events[0]))
    else:
        mave_hgvs = "{}.[{}]".format(
            prefix, ";".join([format_variant(e) for e in events])
        )

    match = single_variant_re.fullmatch(mave_hgvs) or multi_variant_re.fullmatch(
        mave_hgvs
    )
    if not match:
        raise exceptions.HGVSMatchError(
            "Could not validate parsed variant '{hgvs}'.".format(hgvs=mave_hgvs)
        )

    return mave_hgvs


def non_hgvs_columns(columns):
    """
    Takes an iterable of column names and returns a pandas Index object
    containing all entries that are not equal to the HGVS column names.
    The order of the elements is preserved.
    """
    data_columns = [
        x
        for x in columns
        if x != constants.nt_variant_col and x != constants.pro_variant_col
    ]
    return pd.Index(data_columns)


def hgvs_columns(columns):
    """
    Takes an iterable of column names and returns a pandas Index object
    containing all entries that are equal to the HGVS column names.
    The order of the elements is preserved.
    """
    data_columns = [
        x
        for x in columns
        if x == constants.nt_variant_col or x == constants.pro_variant_col
    ]
    return pd.Index(data_columns)