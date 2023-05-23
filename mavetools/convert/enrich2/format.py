import logging

import pandas as pd
from pandas.testing import assert_index_equal

from . import LOGGER, constants, filters, utilities, validators


__all__ = [
    "apply_offset",
    "drop_null",
    "flatten_column_names",
    "get_count_dataframe_by_condition",
    "get_replicate_score_dataframes",
]


logger = logging.getLogger(LOGGER)


def apply_offset(variant, offset, enrich2=None):
    """
    Applies offset to the base position of a HGVS point mutation by
    subtraction. If `enrich2` is not None, then additional validation
    against a wild-type NT and Protein sequence are performed after applicaiton
    of the offset.
    """
    variants = []
    for v in variant.split(","):
        nt_instance = None
        if len(v.strip().split(" ")) == 2:
            nt, pro = v.strip().split(" ")
        elif v.strip()[0] == "p":
            nt, pro = None, v.strip()
        else:
            nt, pro = v.strip(), None

        if nt is not None:
            nt = utilities.NucleotideSubstitutionEvent(nt)
            nt_instance = nt
            nt.position -= offset
            if nt.position < 1:
                raise ValueError("Position after offset {} " "applied to {} is negative.".format(offset, nt.variant))
            if enrich2:
                enrich2.validate_against_wt_sequence(nt.format)
            nt = nt.format

        if pro is not None and pro != "p.=":
            use_brackets = False
            if pro.startswith("(") and pro.endswith(")"):
                pro = pro[1:-1]
                use_brackets = True

            pro = utilities.ProteinSubstitutionEvent(pro)
            if nt_instance is not None:
                pro.position = nt_instance.codon_position()
            else:
                pro_offset = (1, -1)[offset < 0] * (abs(offset) // 3)
                pro.position -= pro_offset

            if enrich2:
                enrich2.validate_against_protein_sequence(pro.format)
            pro = pro.format
            if use_brackets:
                pro = "({})".format(pro)

        variants.append("{} {}".format("" if nt is None else nt, "" if pro is None else pro).strip())

    return ", ".join(variants)


def drop_null(scores_df, counts_df=None):
    """
    Drops null rows and columns. If `counts_df` is not None, then they
    must contain the same index and same variants under the HGVS columns
    `hgvs_nt` and `hgvs_pro`.
    Modification is inplace when `counts_df` is `None`
    Parameters
    ----------
    scores_df : `pd.DataFrame`
        Scores dataframe the columns `hgvs_nt`, `hgvs_pro` and `score`
    counts_df : `pd.DataFrame`
        Counts dataframe containing the columns `hgvs_nt` and `hgvs_pro`
    Raises
    ------
    AssertionError
    ValueError
    Returns
    -------
    tuple[`pd.DataFrame`]
    """
    if counts_df is not None:
        # Apply filters to drop NaN columns and rows before validation
        # Join will discard rows in `counts_df` that have an index that
        # does not appear in `scores_df`. This shouldn't happen since we
        # validate that both indexes are the same.
        assert_index_equal(scores_df.index, counts_df.index)
        validators.validate_datasets_define_same_variants(scores_df, counts_df)
        joint_df = pd.concat(
            objs=[scores_df, counts_df[utilities.non_hgvs_columns(counts_df.columns)]],
            axis=1,
            sort=False,
        )
        assert len(joint_df) == len(counts_df)
        filters.drop_na_columns(joint_df)
        filters.drop_na_rows(joint_df)

        score_columns = list(utilities.hgvs_columns(joint_df.columns)) + list(
            utilities.non_hgvs_columns(scores_df.columns)
        )
        score_columns = [c for c in score_columns if c in joint_df]
        scores_df = joint_df[score_columns]

        count_columns = list(utilities.hgvs_columns(joint_df.columns)) + list(
            utilities.non_hgvs_columns(counts_df.columns)
        )
        count_columns = [c for c in count_columns if c in joint_df]
        counts_df = joint_df[count_columns]

        assert_index_equal(scores_df.index, counts_df.index)
    else:
        filters.drop_na_columns(scores_df)
        filters.drop_na_rows(scores_df)

    return scores_df, counts_df


def flatten_column_names(columns, ordering):
    """
    Takes a column MultiIndex and joins each entry into a single underscore
    separated string based on the indices in the iterable `ordering`.
    Spaces in the MultiIndex entry are replaced with underscores.
    """
    values = [[x[y] for y in ordering] for x in columns]
    cnames = ["_".join(x) for x in values]
    cnames = [x.replace(" ", "_") for x in cnames]
    return cnames


def get_replicate_score_dataframes(store, element=constants.variants_table):
    """
    Return a dictionary of DataFrames, one for each condition in store.
    Store is an open Enrich2 HDF5 file from an Experiment.
    Dictionary keys are condition names.
    """
    condition_dfs = dict()
    idx = pd.IndexSlice

    scores_key = "/main/{}/scores".format(element)
    shared_key = "/main/{}/scores_shared".format(element)
    if scores_key not in store:
        logger.warning("Store is missing key {}. Skipping score file output.".format(scores_key))
        return condition_dfs
    if shared_key not in store:
        logger.warning("Store is missing key {}. Skipping score file output.".format(shared_key))
        return condition_dfs

    for cnd in store["/main/{}/scores".format(element)].columns.levels[0]:
        assert_index_equal(
            store["/main/{}/scores".format(element)][cnd].index,
            store["/main/{}/scores_shared".format(element)][cnd].index,
        )
        condition_dfs[cnd] = store["/main/{}/scores".format(element)].loc[:, idx[cnd, :, :]]
        condition_dfs[cnd].columns = condition_dfs[cnd].columns.levels[1]

        rep_scores = store["/main/{}/scores_shared".format(element)].loc[:, idx[cnd, :, :]]
        rep_scores.columns = flatten_column_names(rep_scores.columns, (2, 1))

        condition_dfs[cnd] = pd.merge(
            condition_dfs[cnd],
            rep_scores,
            how="inner",
            left_index=True,
            right_index=True,
        )
    return condition_dfs


def get_count_dataframe_by_condition(store, cnd, element=constants.variants_table, filtered=None):
    """
    Return a DataFrame corresponding the condition cnd.
    Store is an open Enrich2 HDF5 file from an Experiment.
    filtered is a pandas Index containing variants to include. If it is none,
    the index of the DataFrame's score table for the element is used.
    """
    idx = pd.IndexSlice

    count_key = "/main/{}/counts".format(element)
    if count_key not in store:
        logger.warning("Store is missing key {}. Skipping count file output.".format(count_key))
        return None

    if filtered is None:
        scores_key = "/main/{}/scores".format(element)
        if scores_key not in store:
            logger.warning("Store is missing key {}. Skipping count file output.".format(scores_key))
            return None
        filtered = store["/main/{}/scores".format(element)].index

    # TODO: revisit tests to see if preserving the all-NA rows makes sense
    store_df = store[count_key]
    store_df = store_df.reindex(filtered)
    df = store_df.loc[filtered, idx[cnd, :, :]]
    df.columns = flatten_column_names(df.columns, (1, 2))
    return df
