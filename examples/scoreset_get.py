import attr, os
from pprint import PrettyPrinter

from mavetools.client.client import Client
from mavetools.models.scoreset import ScoreSet


pp = PrettyPrinter(indent=2)
base_url = os.getenv('MAVEDB_BASE_URL', '')
scoreset_urn = 'urn:mavedb:00000001-a-1'

# Generate a new auth_token in your profile and post it here
auth_token = 'AseyaNLLhqv9jAm0joMkq2oqB0bw3GKxTclkT2NtG340RF6CfdM2UC3j8Fv4RpbQ'
client = Client(base_url, auth_token=auth_token) if base_url else Client(auth_token=auth_token)

# GET
scoreset = client.get_model_instance(ScoreSet, scoreset_urn)
pp.pprint(attr.asdict(scoreset))
