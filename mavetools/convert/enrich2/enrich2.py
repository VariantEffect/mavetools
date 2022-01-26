import logging
from abc import ABCMeta, abstractmethod


class Enrich2(metaclass=ABCMeta):
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
        Set as `True` if the input variants contain coding HGVS syntax
        (c or p). Used only in Enrich2.
    input_type : str, optional.
        The MaveDB file type. Can be either 'scores' or 'counts'.
    """
    LOG_MSG = "Writing {elem} {df_type} for condition '{cnd}' to '{path}'."

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

        if is_coding and not abs(offset) % 3 == 0:
            raise ValueError(
                "Enrich2 offset for a coding " "dataset must be a multiple of 3."
            )

    def convert(self):
        """

        Returns
        -------

        """
        if self.input_is_h5:
            input_file = self.load_input_file()
            result = self.parse_input(input_file)
            input_file.close()
            return result
        else:
            return self.parse_tsv_input(self.load_input_file())