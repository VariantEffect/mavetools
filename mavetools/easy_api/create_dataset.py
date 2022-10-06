import pandas as pd
from mavetools.client.client import Client


def create_scoreset(new_scoreset, auth_token, scores_file_path, counts_file_path=None):
    """

    """
    client = Client(auth_token=auth_token)
    scores = pd.read_csv(scores_file_path)
    if counts_file_path is not None:
        counts = pd.read_csv(counts_file_path)
    scoreset_urn = client.create_scoreset(new_scoreset, scores_df=scores, counts_df=counts)
    return scoreset_urn


def create_experiment(new_experiment, auth_token):
    """

    """
    client = Client(auth_token=auth_token)
    experiment_urn = client.create_experiment(new_experiment)
    return experiment_urn
