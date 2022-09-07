from mavetools.client.client import Client


def get_scoreset(scoreset_urn, base_url):
    """

    """
    client = Client(base_url)
    scoreset = client.get_scoreset(scoreset_urn)
    return scoreset


def get_experiment(experiment_urn, base_url):
    """

    """
    client = Client(base_url)
    experiment = client.get_experiment(experiment_urn)
    return experiment

