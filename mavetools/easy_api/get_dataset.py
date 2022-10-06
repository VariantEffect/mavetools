from mavetools.client.client import Client


def get_scoreset(scoreset_urn):
    """

    """
    client = Client()
    scoreset = client.get_scoreset(scoreset_urn)
    return scoreset


def get_experiment(experiment_urn):
    """

    """
    client = Client()
    experiment = client.get_experiment(experiment_urn)
    return experiment

