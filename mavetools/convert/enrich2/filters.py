import logging

import numpy as np

from . import LOGGER, constants, utilities


logger = logging.getLogger(LOGGER)


def drop_na_columns(df):
    """
    Drop columns where all entries are null. Operation is performed in place.

    Parameters
    __________
    df :


    Returns
    _______
    df :
    """
    has_nt_col = constants.nt_variant_col in df.columns
    has_pro_col = constants.pro_variant_col in df.columns

    if has_nt_col:
        nt_all_null = np.all(
            df.loc[:, constants.nt_variant_col].apply(utilities.is_null)
        )
        if nt_all_null:
            df.drop(columns=[constants.nt_variant_col], inplace=True)
    if has_pro_col:
        pro_all_null = np.all(
            df.loc[:, constants.pro_variant_col].apply(utilities.is_null)
        )
        if pro_all_null:
            df.drop(columns=[constants.pro_variant_col], inplace=True)

    # Drop data columns that are all null.
    to_drop = list()
    for cname in utilities.non_hgvs_columns(df.columns):
        if np.all(df.loc[:, cname].isnull()):
            logger.warning(
                "Dropping column '{}' because it contains all null "
                "values".format(cname)
            )
            to_drop.append(cname)
    if len(to_drop) > 0:
        df.drop(columns=to_drop, inplace=True)

    return df


def drop_na_rows(df):
    """
    Drop rows where all non-HGVS entries are null. Operation is performed in
    place.

    Parameters
    __________
    df

    Returns
    _______
    df
    """
    null_rows = df.loc[:, utilities.non_hgvs_columns(df.columns)].isnull().all(axis=1)
    if sum(null_rows) > 0:
        logger.warning(
            "Dropping {} rows that contain all null values".format(sum(null_rows))
        )
        df.drop(index=df.index[null_rows], inplace=True)

    return df