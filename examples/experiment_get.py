import attr, os
from pprint import PrettyPrinter

from mavetools.client.client import Client
from mavetools.models.experiment import Experiment


pp = PrettyPrinter(indent=2)
base_url = os.getenv('MAVEDB_BASE_URL', '')
experiment_urn = 'urn:mavedb:00000001-a'

# Generate a new auth_token in your profile and post it here
auth_token = 'AseyaNLLhqv9jAm0joMkq2oqB0bw3GKxTclkT2NtG340RF6CfdM2UC3j8Fv4RpbQ'
client = Client(base_url, auth_token=auth_token) if base_url else Client(auth_token=auth_token)


experiment = client.get_model_instance(Experiment, experiment_urn)
pp.pprint(attr.asdict(experiment))
