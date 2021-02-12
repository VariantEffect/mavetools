import attr, os
from pprint import PrettyPrinter

from mavetools.client.client import Client
from mavetools.models.licence import Licence
from mavetools.models.scoreset import NewScoreSet, NewScoreSetRequest, ScoreSet
from mavetools.models.target import NewTarget, ReferenceGenome, ReferenceMap, SequenceOffset


pp = PrettyPrinter(indent=2)
base_url = os.getenv('MAVEDB_BASE_URL', '')
scoreset_urn = 'urn:mavedb:00000001-a-1'
experiment_urn = 'urn:mavedb:00000001-a'
# Generate a new auth_token in your profile and post it here
auth_token = 'jErFdtOCBjxmiYllKwEtzTqifAo2PSG9yaVJux8SJcPEMxCFkFbv5wQQkuB7nNUv'
client = Client(base_url, auth_token=auth_token) if base_url else Client(auth_token=auth_token)


# GET
scoreset = client.get_model_instance(ScoreSet, scoreset_urn)
pp.pprint(attr.asdict(scoreset))



# POST
# Change this dir string as needed. It's currently configured for running
# inside a Docker container that mounts the home directory as a volume.
test_file_dir = '/mavetools/tests/test_upload_scoreset/test_files'
new_scoreset = NewScoreSet(
    title='test_title',
    short_description='test_short_description',
    abstract_text='test_abstract_text',

    experiment=experiment_urn,
    score_data=f"{test_file_dir}/test_score_data.csv",
    count_data=f"{test_file_dir}/test_count.csv",
    meta_data=f"{test_file_dir}/test_metadata.json",
    licence=Licence(short_name='CC BY 4.0'),
)
new_scoreset_request = NewScoreSetRequest(
    scoreset=new_scoreset,
    target=NewTarget(
        name='test_target_name',
        type='Protein coding',
        sequence_type='Infer',
        fasta_file=f"{test_file_dir}/test_fasta_file.fasta"
    ),
    uniprot=SequenceOffset(offset=1, identifier='P63165'),
    ensembl=SequenceOffset(offset=1, identifier='ENSG00000116030'),
    refseq=SequenceOffset(offset=1, identifier='NM_001005781.1'),
    reference_maps=[
        ReferenceMap(genome=ReferenceGenome(short_name='hg16'))
    ]
)
client.post_model_instance(new_scoreset_request)
