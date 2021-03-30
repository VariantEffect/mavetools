import pandas as pd
from csv import reader


def df_to_pandas(df):
    """
    This function converts a dataframe downloaded from MaveDB into a pandas dataframe and returns
    the new data frame. The user has the option to drop the accession numbers.

    Parameters
    ----------
    df (string): the filename.csv to be converted to pandas df

    Returns
    -------
    df_pandas (pandas.df): the converted df

    Raises
    ______
    TypeError:
        if argument is not string
        if no arguments are passed
    ValueError:
        if first argument is not in filename.csv format
    """
    # check for TypeError
    # if first argument is not string
    #if not isinstance(df, str):
    #    raise TypeError("df must be string")
    # if no arguments are passed
    #if df is None:  #************** remove this
    #    raise TypeError("no arguments passed")

    # check for ValueError
    # if df is not in filename.csv format
    #if not df.endswith(".csv"):
    #    raise ValueError("df must be csv file")

    # convert df to pandas df
    df_pandas = pd.read_csv(df, skiprows=4)

    # return df_pandas and meta_dict
    return df_pandas


def get_meta_data(df):
    """
    This function takes in a df and returns a dictionary with the formatted metadata.

    Parameters
    ----------
    df (string): the filename.csv to be converted to pandas df

    Returns
    -------
    meta_dict (dictionary): dictionary containing formatted metadata (comments) at the top of df

    Raises
    ______

    """

    # store metadata in dictionary
    # declare dictionary
    meta_dict = {}

    # open csv file in read mode
    with open("urn mavedb 00000001-a-1_scores.csv", "r") as mave_data:  # *********** dNOT hardcode this
        csv_reader = reader(mave_data)
        # read line by line and check for metadata
        for row in csv_reader:
            if not row[0].startswith("#"):  # we have found all the metadata
                break
            else:
                # get the metadata at that index
                get_meta = row[0]
                # remove '#' symbol
                get_meta = get_meta[2:]
                # break string into key/value pairs
                split_meta = get_meta.split(":", 1)
                # add to meta_dict
                meta_dict.update({split_meta[0]: split_meta[1][1:]})

    # return df_pandas and meta_dict
    return meta_dict
