import _io
import pandas as pd
from csv import reader


def df_to_pandas(df, ret_meta=False):
    """
    This function converts a dataframe downloaded from MaveDB into a pandas dataframe and returns
    the new data frame and the metadata formatted as a dictionary if desired. The user has the
    option to drop the accession numbers.

    Parameters
    ----------
    df : string
        the filename.csv to be converted to pandas df
    ret_meta : bool
        True if user wants metadata,
        default value = False

    Returns
    -------
    df_pandas : pandas.df
        the converted df
    meta_dict : dict
        dictionary containing formatted metadata (comments) at the top of df

    Raises
    ______
    TypeError
        if df is not string
    TypeError
        if ret_meta is not bool
    ValueError
        if no arguments are passed
    ValueError
        if first argument is not in filename.csv format
    """
    # check for TypeError
    # if df is not string
    if not isinstance(df, str) and not isinstance(
        df, _io.StringIO
    ):  # account for StringIO object in test file
        raise TypeError("df must be string in filename.csv format")
    # if ret_meta is not bool
    if not isinstance(ret_meta, bool):
        raise TypeError("ret_meta must be bool")

    # check for ValueError
    # if df is not in filename.csv format
    if not isinstance(df, _io.StringIO) and not df.endswith(
        ".csv"
    ):  # account for StringIO object in test file
        raise ValueError("df must be csv file")

    # store metadata in dictionary
    # declare dictionary
    meta_dict = {}

    # open csv file in read mode, if it's a csv file
    if isinstance(df, str):  # df is a filename.csv file
        mave_data = open(df, "r")
    elif isinstance(df, _io.StringIO):  # df is a StringIO object for testing
        mave_data = df

    # get metadata
    for row in range(4):
        get_meta = mave_data.readline()
        # remove '#' symbol
        get_meta = get_meta[2:]
        # break string into key/value pairs
        split_meta = get_meta.split(":", 1)
        # add to meta_dict
        meta_dict.update({split_meta[0]: split_meta[1][1:]})

    # convert df to pandas df
    df_pandas = pd.read_csv(df)

    # close file
    mave_data.close()

    # return df_pandas and meta_dict
    if not ret_meta:
        return df_pandas
    else:
        return df_pandas, meta_dict
