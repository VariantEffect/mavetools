import os
import re
import logging
import numpy as np
from abc import ABCMeta, abstractmethod

# from hgvsp import is_multi
from mavehgvs import Variant
from fqfa.constants.iupac.protein import AA_CODES
from fqfa.validator.validator import dna_bases_validator

from . import LOGGER, utilities, constants

logger = logging.getLogger(LOGGER)

__all__ = ["BaseProgram"]


class BaseProgram(metaclass=ABCMeta):
    """
    Convert an input file to MaveDB_ compliant counts or scores files.
    Attributes
    ----------
    src : str
        Source file to convert.
    wt_sequence : str
        A DNA wild-type sequence used for validation and inference of
        missing variants.
    offset : int, optional.
        The number of bases to clip in `wt_seq`. The resulting sequence
        should be the coding sequence analyzed in your input file.
    dst : str, optional.
        Directory to save the output to. Inferred as input directory if not
        specified.
    one_based : bool, optional.
        Set to `True` if input positions are `1-based` relative to the
        the wild-type sequence specified.
    skip_header_rows : int
        The number of lines to skip at the start of the input file. Only
        applicable to `Excel` and `TSV` files.
    skip_header_rows : int
        The number of lines to skip at the end of the input file. Only
        applicable to `Excel` and `TSV` files.
    sheet_name : str, optional.
        Name of the sheet to convert in an excel file.
    score_column : str, optional.
        The name of the column in the input column which you would like to
        set as the MaveDB score column. Ignored if `input_type` is 'counts'.
        Used when input is tsv.

    hgvs_column : str, optional.
        The name of the column in the input column which you would like to
        set as the MaveDB hgvs_nt column. Used when input is an Enrich2 tsv
        dataset.

    is_coding : bool, optional.
        Set as `True` if the input variants are contain coding HGVS syntax
        (c or p). Used only in Enrich2.
    input_type : str, optional.
        The MaveDB file type. Can be either 'scores' or 'counts'.
    """

    def __init__(
        self,
        src,
        wt_sequence,
        offset=0,
        dst=None,
        one_based=True,
        skip_header_rows=0,
        skip_footer_rows=0,
        score_column="score",
        hgvs_column="hgvs",
        input_type=None,
        sheet_name=None,
        is_coding=True,
    ):
        # Check the input is a readable file.
        self.src = os.path.normpath(os.path.expanduser(src))
        logger.info("Checking read permission for '{}'".format(self.src))
        os.access(self.src, os.R_OK)

        src_filename, ext = os.path.splitext(os.path.split(src)[1])
        self.src_filename = src_filename
        self.dst_filename = "mavedb_{}.csv".format(re.sub(r"\s+", "_", src_filename))
        self.ext = ext.lower()

        # Set directory as the same directory as the input file if not provided
        # For HDF5 files, create a further directory with the same name
        # as the input file since there will be multiple output files.
        self.dst = dst
        if self.dst is None:
            dst, _ = os.path.split(src)
            if self.ext.lower() == ".h5":
                dst = os.path.normpath(os.path.join(os.path.expanduser(dst), self.src_filename))
            self.dst = dst
        else:
            self.dst = os.path.normpath(os.path.expanduser(dst))
        # Create directory tree if it does not exist and check for
        # read and write permissions.
        if not os.path.isdir(self.dst):
            logger.info("Creating directory '{}'".format(self.dst))
            os.makedirs(self.dst, exist_ok=True)
        logger.info("Checking write permission to directory '{}'".format(self.dst))
        os.access(self.dst, os.W_OK)

        self.is_coding = is_coding
        self.skip_header_rows = skip_header_rows
        self.skip_footer_rows = skip_footer_rows
        self.sheet_name = sheet_name
        self.score_column = score_column
        self.hgvs_column = hgvs_column
        self.input_type = input_type
        self.one_based = one_based

        # Initialize sequence information.
        self._wt_sequence = None
        self.codons = None
        self.protein_sequence = None

        self.offset = offset
        self.wt_sequence = wt_sequence

    @property
    def wt_sequence(self):
        return self._wt_sequence

    @wt_sequence.setter
    def wt_sequence(self, seq):
        seq = str(seq).upper()
        # Initialize sequence information.
        if dna_bases_validator(seq) is None:
            raise ValueError("{} is not a valid DNA sequence.".format(seq))
        if self.is_coding:
            self.protein_sequence = utilities.translate_dna(seq, offset=0)
            self.codons = list(utilities.slicer(seq, 3))
        self._wt_sequence = seq

    @property
    def extension(self):
        return self.ext.lower()

    @property
    def input_is_h5(self):
        return self.ext.lower() == ".h5"

    @property
    def input_is_tsv(self):
        return self.ext.lower() == ".tsv"

    @property
    def input_is_scores_based(self):
        return self.input_type == constants.score_type

    @property
    def output_directory(self):
        return os.path.normpath(os.path.expanduser(self.dst))

    @property
    def output_file(self):
        return os.path.normpath(os.path.join(self.output_directory, self.dst_filename))

    def convert(self):
        """
        Runs `parse_input` and saves the Mavedb-compliant result to file.
        """
        logger.info("Processing file {}".format(self.src))
        mave_df = self.parse_input(self.load_input_file())
        logger.info("Writing to {}".format(self.output_file))
        mave_df.to_csv(self.output_file, sep=",", index=None, na_rep=np.NaN)

    @abstractmethod
    def load_input_file(self):
        pass  # pragma: no cover

    @abstractmethod
    def parse_input(self, df):
        pass  # pragma: no cover

    @abstractmethod
    def parse_row(self, row):
        pass  # pragma: no cover

    def validate_against_wt_sequence(self, variant):
        # TODO
        # this validation code exists elsewhere in mavetools
        """
        Checks that the reference base in a substitution variant matches that
        in the wild-type sequence provided.
        Parameters
        ----------
        variant : str
            A nucleotide substitution variant with valid HGVS_ syntax.
        """
        # if is_multi(variant):
        if Variant(variant).is_multi_variant():
            _ = [self.validate_against_wt_sequence(v) for v in utilities.split_variant(variant)]
            return

        if variant in constants.special_variants:
            return

        variant = utilities.NucleotideSubstitutionEvent(variant)
        if variant.silent:
            return

        zero_based_pos = variant.position - int(self.one_based)
        if zero_based_pos < 0:
            raise IndexError(
                ("Encountered a negative position in {}. " "Positions might not be one-based.").format(
                    variant, self.one_based
                )
            )

        if zero_based_pos >= len(self.wt_sequence):
            raise IndexError(
                "Position {} (index {}) extends beyond the maximum index {} in "
                "the wild-type sequence {} with length {}.".format(
                    zero_based_pos + int(self.one_based),
                    zero_based_pos,
                    variant,
                    len(self.wt_sequence) - 1,
                    self.wt_sequence,
                    len(self.wt_sequence),
                )
            )

        wt_ref_nt = self.wt_sequence[zero_based_pos]
        if variant.ref != wt_ref_nt:
            raise ValueError(
                "Reference base '{base}' at 1-based position {pos} in the "
                "wild-type sequence does not match the reference base '{ref}' "
                "suggested in variant '{variant}'.".format(
                    pos=zero_based_pos + 1,
                    base=wt_ref_nt,
                    variant=variant,
                    ref=variant.ref,
                )
            )

    def validate_against_protein_sequence(self, variant):
        """
        Checks that the reference amino acid in a substitution variant matches
        that in the translated wild-type sequence provided.
        Parameters
        ----------
        variant : str
            A protein substitution variant with valid HGVS_ syntax.
        """
        # if is_multi(variant):
        if Variant(variant).is_multi_variant():
            _ = [self.validate_against_protein_sequence(v) for v in utilities.split_variant(variant)]
            return

        if variant in constants.special_variants or "p.=" in variant:
            return

        variant = utilities.ProteinSubstitutionEvent(variant)

        if variant.position > len(self.protein_sequence):
            raise IndexError(
                "Position {} in {} "
                "extends beyond the maximum index {} in the translated "
                "wild-type sequence {} with length {}.".format(
                    variant.position,
                    variant,
                    len(self.protein_sequence) - 1,
                    self.protein_sequence,
                    len(self.protein_sequence),
                )
            )

        wt_aa = AA_CODES[self.protein_sequence[variant.position - 1]]
        if variant.ref != wt_aa:
            raise ValueError(
                "Reference AA '{aa}' at 1-based position {pos} in the "
                "translated protein sequence {seq} does not match the "
                "reference AA '{ref}' suggested in variant '{variant}'.".format(
                    pos=variant.position,
                    aa=wt_aa,
                    variant=variant,
                    ref=variant.ref,
                    seq=self.protein_sequence,
                )
            )
