import pandas as pd
import asyncio
from aiohttp import ClientSession, TCPConnector

from mavetools.client.client import Client
from mavetools.client_async.client import Client as ClientAsync


def create_experiment(new_experiment, auth_token):
    """

    """
    client = Client(auth_token=auth_token)
    experiment_urn = client.create_experiment(new_experiment)
    return experiment_urn


async def create_experiments(experiment_list, auth_token):
    """

    """
    client = ClientAsync(auth_token=auth_token)
    async with ClientSession(connector=TCPConnector(ssl=client.sslcontext)) as client.session:
        r = await asyncio.gather(*[client.create_dataset(experiment, "experiments") for experiment in experiment_list])
    return r


def create_scoreset(new_scoreset, auth_token, scores_file_path, counts_file_path=None):
    """

    """
    client = Client(auth_token=auth_token)
    scores = pd.read_csv(scores_file_path)
    if counts_file_path is not None:
        counts = pd.read_csv(counts_file_path)
    scoreset_urn = client.create_scoreset(new_scoreset, scores_df=scores, counts_df=counts)
    return scoreset_urn


async def create_scoresets(scoreset_list, auth_token):
    """

    """
    client = ClientAsync(auth_token=auth_token)
    async with ClientSession(connector=TCPConnector(ssl=client.sslcontext)) as client.session:
        r = await asyncio.gather(*[client.create_dataset(scoreset[0],
                                                         "scoresets",
                                                         scoreset[1],
                                                         scoreset[2]) for scoreset in scoreset_list])
    return r
