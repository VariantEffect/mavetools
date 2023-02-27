import asyncio
from aiohttp import ClientSession, TCPConnector

from mavetools.client.client import Client as Client
from mavetools.client_async.client import Client as ClientAsync


def get_experiment(experiment_urn):
    """

    """
    client = Client()
    experiment = client.get_experiment(experiment_urn)
    return experiment


async def get_experiments(urn_list):
    """

    """
    client = ClientAsync()
    async with ClientSession(connector=TCPConnector(ssl=client.sslcontext)) as client.session:
        r = await asyncio.gather(*[client.get_dataset("experiments", urn) for urn in urn_list])
    return r


def get_scoreset(scoreset_urn):
    """

    """
    client = Client()
    scoreset = client.get_scoreset(scoreset_urn)
    return scoreset


async def get_scoresets(urn_list):
    """

    """
    client = ClientAsync()
    async with ClientSession(connector=TCPConnector(ssl=client.sslcontext)) as client.session:
        r = await asyncio.gather(*[client.get_dataset("scoresets", urn) for urn in urn_list])

    return r
