import pandas as pd
from csv import reader

def df_to_pandas(df, drop_accession=False):
    """
    This function converts a dataframe downloaded from MaveDB into a pandas dataframe.
    The user has the option to drop the accession numbers.

    Parameters
    ----------
    df (string): the filename.csv to be converted to pandas df
    drop_accession (bool): True if user wants to drop the accession numbers
        default value = False

    Returns
    -------
    df_pandas (pandas.df): the converted df
    meta_dict (dictionary): dictionary containing formatted metadata (comments) at the top of df

    Raises
    ______
    ValueError: if first argument is not in filename.csv format
    TypeError: if second argument is not bool
    TypeError: if no arguments are passed
    """

    # if df is not in filename.csv format
    if not df.endswith(".csv"):
        raise ValueError("df must be csv file")
    # if drop_accession is not bool
    if not isinstance(drop_accession, bool):
        raise TypeError("drop_accession must be boolean value")
    # if no arguments are passed
    if df is None:
        raise TypeError("no arguments passed")

    # convert df to pandas df
    df_pandas = pd.read_csv(df, skiprows=4)

    # remove accession numbers if drop_accession = True
    if drop_accession:
        # make sure accession numbers exist
        if "accession" in df_pandas:
            del df_pandas["accession"]

    # store metadata in dictionary
    # declare dictionary
    meta_dict = {}

    # open csv file in read mode
    with open("urn mavedb 00000001-a-1_scores.csv", "r") as mave_data:
        csv_reader = reader(mave_data)
        # read line by line and check for metadata
        for row in csv_reader:
            if row[0].startswith("#") == False: # we have found all the metadata
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
    return df_pandas, meta_dict














