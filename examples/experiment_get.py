import attr, os
from pprint import PrettyPrinter

from mavetools.client.client import Client
from mavetools.models.experiment import Experiment


pp = PrettyPrinter(indent=2)  # displayes results in readable format
# check environment variables and see if variable named MAVEDB_BASE_URL exists and return value
# if the value does not exist, an empty string is returned instead
base_url = os.getenv("MAVEDB_BASE_URL", "")
experiment_urn = "urn:mavedb:00000001-a"  # the urn of the experiment we want to get

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

# using the client object, GET the model instance of an Experiment with a particular urn
# GET retrieves a resource from the server via the appropriate API endpoint
experiment = client.get_model_instance(Experiment, experiment_urn)
# display results
pp.pprint(attr.asdict(experiment))
