import attr, os
from pprint import PrettyPrinter

from mavetools.client.client import Client
from mavetools.models.scoreset import ScoreSet


pp = PrettyPrinter(indent=2)  # for formatting output
# check environment variables and see if variable named MAVEDB_BASE_URL exists and return value
# if the value does not exist, an empty string is returned instead
base_url = os.getenv("MAVEDB_BASE_URL", "")
scoreset_urn = "urn:mavedb:00000001-a-1"  # the urn of the experiment we want to get

# Generate a new auth_token in your profile and post it here
auth_token = "AseyaNLLhqv9jAm0joMkq2oqB0bw3GKxTclkT2NtG340RF6CfdM2UC3j8Fv4RpbQ"
# auth_token =
# if the base url exists, the client object is instantiated with that value
# otherwise the client object is instantiated with default value which points to localhost
client = (
    Client(base_url, auth_token=auth_token)
    if base_url
    else Client(auth_token=auth_token)
)

# GET
scoreset = client.get_model_instance(ScoreSet, scoreset_urn)
pp.pprint(attr.asdict(scoreset))
