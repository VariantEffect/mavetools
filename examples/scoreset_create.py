import attr, os

from mavetools.client.client import Client
from mavetools.models.licence import Licence
from mavetools.models.scoreset import NewScoreSet, NewScoreSetRequest, ScoreSet
from mavetools.models.target import (
    NewTarget,
    ReferenceGenome,
    ReferenceMap,
    SequenceOffset,
)

# check environment variables and see if variable named MAVEDB_BASE_URL exists and return value
# if the value does not exist, an empty string is returned instead
from mavetools.validators.validate import validate_all

base_url = os.getenv("MAVEDB_BASE_URL", "")
# the urn of the scoreset and the experiment that that scoreset belongs to
# what happens if the urn of the scoreset already exists?
scoreset_urn = "urn:mavedb:00000001-a-1"
experiment_urn = "urn:mavedb:00000001-a"

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

# POST
# Change this dir string as needed. It's currently configured for running
# inside a Docker container that mounts the home directory as a volume.
test_file_dir = "/mavetools/tests/test_upload_scoreset/test_files"

# validate
# validate_all(countfile=f"{test_file_dir}/test_count.csv",
#             scorefile=f"{test_file_dir}/test_score_data.csv",
#             scorejson=None)

new_scoreset = NewScoreSet(
    title="test_title",
    short_description="test_short_description",
    abstract_text="test_abstract_text",
    experiment=experiment_urn,
    score_data=f"{test_file_dir}/test_score_data.csv",
    count_data=f"{test_file_dir}/test_count.csv",
    meta_data=f"{test_file_dir}/test_metadata.json",
    licence=Licence(short_name="CC BY 4.0"),
    sra_ids=["SRP109119"],
    pubmed_ids=["23035249"],
    doi_ids=["10.1038/s41467-019-11526-w"],
)
new_scoreset_request = NewScoreSetRequest(
    scoreset=new_scoreset,
    target=NewTarget(
        name="test_target_name",
        type="Protein coding",
        sequence_type="Infer",
        fasta_file=f"{test_file_dir}/test_fasta_file.fasta",
    ),
    uniprot=SequenceOffset(offset=1, identifier="P63165"),
    ensembl=SequenceOffset(offset=1, identifier="ENSG00000116030"),
    refseq=SequenceOffset(offset=1, identifier="NM_001005781.1"),
    reference_maps=[ReferenceMap(genome=ReferenceGenome(short_name="hg16"))],
)
client.post_model_instance(new_scoreset_request)
