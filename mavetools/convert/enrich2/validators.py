import logging
from abc import ABCMeta, abstractmethod

from mavehgvs.patterns import dna

import numpy as np
import pandas as pd
from numpy.testing import assert_array_equal

from tqdm import tqdm

from joblib import Parallel, delayed

from mavetools.convert.enrich2 import constants, utilities, exceptions, LOGGER


logger = logging.getLogger(LOGGER)


class ValidationBackend(metaclass=ABCMeta):
    """
    Validation backend which provides the interface `validate` for validating
    HGVS_ variants.
    """

    @abstractmethod
    def validate(self, variant):
        pass  # pragma: no cover


class HGVSPatternsBackend(ValidationBackend):
    """
    Backend using the regex based validation in `hgvsp`. Fast but may be
    too strict in specific cases.
    """

    def validate(self, variant):
        """
        Validates a HGVS_ variant using the regex patterns in `hgvsp`.
        Parameters
        ----------
        variant : str
            HGVS formatted variant string.
        Returns
        -------
        str
        """
        if variant in constants.special_variants:
            return variant
        #single_match = hgvsp.single_variant_re.fullmatch(variant)
        single_match = dna.dna_single_variant.fullmatch(variant)
        #multi_match = hgvsp.multi_variant_re.fullmatch(variant)
        multi_match = dna.dna_multi_variant.fullmatch(variant)
        if not (single_match or multi_match):
            raise exceptions.HGVSValidationError(
                "'{}' is not valid HGVS syntax.".format(variant)
            )
        return variant


def validate_variants(
    variants, validation_backend=None, n_jobs=1, verbose=0, backend="multiprocessing"
):
    """
    Validate each variant's HGVS_ syntax.
    Parameters
    ----------
    variants : list[str]
        Variant HGVS_ representations.
    validation_backend : ValidationBackend
        A parsing backend implementing `validate`.
    n_jobs : int, optional
        Number of jobs to run in parallel.
    verbose : int, optional
        Joblib's verbosity level.
    backend : str, optional
        Parallel backend to use. Defaults to `multiprocessing`.
    Returns
    -------
    list[Union[str, SequenceVariant]]
        Formatted and validated variants.
    """
    if validation_backend is None:
        validation_backend = HGVSPatternsBackend()
    return Parallel(n_jobs=n_jobs, verbose=verbose, backend=backend)(
        delayed(validation_backend.validate)(variant) for variant in variants
    )


def validate_has_column(df, column):
    """Validates that a `DataFrame` contains `column` in it's columns."""
    if column not in df.columns:
        raise KeyError(
            "Missing column '{}'. Existing columns are {}.".format(
                column, ", ".join(df.columns)
            )
        )


def validate_columns_are_numeric(df):
    """Checks non-hgvs columns for float or int data."""
    for column in df.columns:
        if column in [constants.pro_variant_col, constants.nt_variant_col]:
            continue
        else:
            if not (
                np.issubdtype(df.dtypes[column], np.floating)
                or np.issubdtype(df.dtypes[column], np.integer)
            ):
                raise TypeError(
                    "Expected only float or int data columns. Got {}.".format(
                        str(df.dtypes[column])
                    )
                )


def validate_hgvs_uniqueness(df: pd.DataFrame, cname: str) -> None:
    """Validate that the HGVS column entries are unique.
    Parameters
    ----------
    df : pd.DataFrame
        The data frame to validate.
    cname : str
        The column name for the HGVS strings.
    Returns
    -------
    None
    Raises
    ------
    KeyError
        If cname is not a column name in df.
    ValueError
        If there are non-unique values in cname (not including None).
    """
    try:
        values = df[cname]
    except KeyError:
        raise KeyError(f"invalid column name '{cname}'")
    else:
        dup_counts = values.value_counts()
        dups = dup_counts[dup_counts > 1].index
        if len(dups) > 0:
            dup_error_string = ", ".join(dups[: constants.MAX_ERROR_VARIANTS])
            if len(dups) > constants.MAX_ERROR_VARIANTS:
                dup_error_string += ", ..."
            raise ValueError(
                f"found {len(dups)} duplicate HGVS strings in '{cname}': {dup_error_string}"
            )


def validate_datasets_define_same_variants(scores_df, counts_df):
    """
    Checks if two `pd.DataFrame` objects parsed from uploaded files
    define the same variants.
    Parameters
    ----------
    scores_df : `pd.DataFrame`
        Scores dataframe parsed from an uploaded scores file.
    counts_df : `pd.DataFrame`
        Scores dataframe parsed from an uploaded counts file.
    """

    scores_columns = [c for c in scores_df.columns if c in constants.variant_columns]
    counts_columns = [c for c in counts_df.columns if c in constants.variant_columns]
    if scores_columns != counts_columns:
        raise AssertionError(
            "Dataframes define different hgvs columns. "
            "Scores defines '{}' and counts defines '{}'.".format(
                ", ".join(scores_columns), ", ".join(counts_columns)
            )
        )

    if constants.nt_variant_col in scores_columns:
        scores_nt = scores_df[constants.nt_variant_col].values
        counts_nt = counts_df[constants.nt_variant_col].values
        try:
            assert_array_equal(scores_nt, counts_nt)  # Treats np.NaN as equal
        except AssertionError:
            not_equal_selector = scores_nt != counts_nt
            neq_list = [
                "{} ({})".format(x, y)
                for x, y in zip(
                    scores_nt[not_equal_selector], counts_nt[not_equal_selector]
                )
                if (x is not np.NaN) and (y is not np.NaN)
            ]
            raise AssertionError(
                "Scores and counts do not define the same "
                "nucleotide variants: {}.".format(", ".join(neq_list))
            )

    if constants.pro_variant_col in scores_columns:
        scores_pro = scores_df[constants.pro_variant_col].values
        counts_pro = counts_df[constants.pro_variant_col].values
        try:
            assert_array_equal(scores_pro, counts_pro)  # Treats np.NaN as equal
        except AssertionError:
            not_equal_selector = scores_pro != counts_pro
            neq_list = [
                "{} ({})".format(x, y)
                for x, y in zip(
                    scores_pro[not_equal_selector], counts_pro[not_equal_selector]
                )
                if (x is not np.NaN) and (y is not np.NaN)
            ]
            raise AssertionError(
                "Scores and counts do not define the same protein variants: "
                "{}.".format(", ".join(neq_list))
            )


def validate_mavedb_compliance(df, df_type):
    """Runs MaveDB compliance checks."""
    tqdm.pandas(desc="Validating variants")

    has_nt_col = constants.nt_variant_col in df.columns
    has_pro_col = constants.pro_variant_col in df.columns
    if not has_nt_col and not has_pro_col:
        raise ValueError(
            "Dataframe must define either '{}', '{}' or both.".format(
                constants.nt_variant_col, constants.pro_variant_col
            )
        )

    primary_col = None
    if has_nt_col:
        defines_nt = not all(
            df.loc[:, constants.nt_variant_col].progress_apply(utilities.is_null)
        )
        if defines_nt:
            primary_col = constants.nt_variant_col

    if has_pro_col and primary_col is None:
        defines_pro = not all(
            df.loc[:, constants.pro_variant_col].progress_apply(utilities.is_null)
        )
        if defines_pro:
            primary_col = constants.pro_variant_col

    if primary_col is None:
        raise ValueError(
            "Neither '{}' or '{}' defined any variants.".format(
                constants.nt_variant_col, constants.pro_variant_col
            )
        )

    null_primary = df.loc[:, primary_col].progress_apply(utilities.is_null)
    if any(null_primary):
        raise ValueError(
            "Primary column (inferred as '{}') cannot "
            "contain the null values {} (case-insensitive).".format(
                primary_col, "NaN, Na, None, whitespace, Undefined"
            )
        )

    try:
        validate_hgvs_uniqueness(df, primary_col)
    except ValueError as e:
        # allow duplicates for protein primary
        # convert error to warning
        if primary_col == constants.pro_variant_col:
            logger.warning(e)
        else:
            raise e

    validate_columns_are_numeric(df)
    if df_type == constants.score_type:
        validate_has_column(df, "score")
    return df