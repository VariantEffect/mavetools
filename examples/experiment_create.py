import os

from mavetools.client.client import Client
from mavetools.models.experiment import NewExperiment

base_url = os.getenv('MAVEDB_BASE_URL', '')
# experimentset_urn = 'tmp:jCICvwLCntIuKIsf'

# Generate a new auth_token in your profile and post it here
auth_token = 'AseyaNLLhqv9jAm0joMkq2oqB0bw3GKxTclkT2NtG340RF6CfdM2UC3j8Fv4RpbQ'
client = Client(base_url, auth_token=auth_token) if base_url else Client(auth_token=auth_token)

# POST
new_experiment = NewExperiment(
    # experimentset=experimentset_urn,
    title='exp_test_title',
    short_description='exp_test_short_description',
    abstract_text='test_abstract_text',
    method_text='test_method_text',

    sra_ids=['SRP109119'],
    pubmed_ids=['23035249'],
    doi_ids=['10.1038/s41467-019-11526-w'],
)
client.post_model_instance(new_experiment)
