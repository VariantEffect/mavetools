POST scoreset
=============

create a scoreset in MaveDB via the API
---------------------------------------

To begin, import the modeules below.

.. code:: ipython3

    import attr, os
    from mavetools.client.client import Client
    from mavetools.models.licence import Licence
    from mavetools.models.scoreset import NewScoreSet, NewScoreSetRequest, ScoreSet
    from mavetools.models.target import NewTarget, ReferenceGenome, ReferenceMap, SequenceOffset
    from mavetools.validators.validate import validate_all

Here your base_url is set to localhost, http://127.0.0.1:8000/api/. This
default funcionality is what you would want to use when working with a
local instance of MaveDB (e.g., a development branch). If working with
production mavedb you would set base url to https://www.mavedb.org/api/.

In the cell below, comment out the base_url you will not be using.

.. code:: ipython3

    base_url = 'http://127.0.0.1:8000/api/'
    #base_url = 'https://www.mavedb.org/api/'

Set the name value of experiment_urn to the urn of the experiment where
the scoreset belongs.

.. code:: ipython3

    experiment_urn = 'tmp:WRe9wTCdGxKKQV4a'

Next, you will need an auth_token to make POST requests to MaveDB. If
you have one, substitute it in the example provided below. If you need
one, please follow these instructions:

::

   1. go to https://www.mavedb.org
   2. login using your ORCID ID
   3. go to settings
   4. generate new auth token
   5. copy auth token and pase it in the auth_token field below

.. code:: ipython3

    # Generate a new auth_token in your profile and post it here
    auth_token = 'R2skRbpBD3Rsf5dNHoQxDZevdEE74T5lCKMFyBhBwwPFH4ZfTrxDz7TZ0kbFLtEZ'

Here you instantiate the Client object. The Client object is the object
by which the POST request is performed. The client object is
instantiated with the value of base_url provided earlier, so make sure
that is up-to-date. If base_url does not exist, base_url is defaulted to
localhost, http://127.0.0.1:8000/api/.

.. code:: ipython3

    client = Client(base_url, auth_token=auth_token) if base_url else Client(auth_token=auth_token)

test_file_dir is the path to the directory in which the files needed for
making a scoreset POST resquest exist. The required files are as
follows:

1. abstract.md
2. method.md
3. test_count.csv
4. test_fasta_file.fasta
5. test_metadata.json
6. test_score_data.csv

For the abstrct.md and method.md files, simply paste your content into
these files.

The test\_ files above will be replaced by your own files. You can do
this in two ways, replace the name of your files to correspond with the
above files and replace the files in the directory listed below
(recommended). Or, you can put your files with their current name in the
directory, just ensure that you change the names accoringly when
instantiating the NewScoreSet and NewScoreSetRequest onjects later in
this module.

We have an example directory in mavetools that holds the files of
interest. Though this directory can exist anywhere on your computer, you
must put the correct path to that directory as the value to
test_file_dir.

.. code:: ipython3

    # here is an example if your copy of mavetools exists within Pycharm
    test_file_dir = '~/PycharmProjects/mavetools/tests/test_upload_scoreset/test_files'

Here we want to valide the data we want to POST, run the validation code
here.

If you get an error message here, DO NOT UPLOAD, as your upload will
fail serverside as well. Instead resolve the error, run script again,
confirm error has been resolved, then upload.

.. code:: ipython3

    #validate
    #validate_all(count_data=f"{test_file_dir}/test_count.csv", 
    #             score_data=f"{test_file_dir}/test_score_data.csv", 
    #             scorejson=None)

Instantiate the NewScoreSet object and assign it to the new_scoreset
veriable. You must substitute the attribute values for your scoreset’s
attribute values.

.. code:: ipython3

    #with open(f"{test_file_dir}/abstract.md") as handle:
    #    test_abstract_text = handle.read()
    
    # substitute each attribute for your scoreset attributes
    new_scoreset = NewScoreSet(
        title='test_title',
        short_description='test_short_description',
        #abstract_text=test_abstract_text,
        abstract_text="abstract",
    
        experiment=experiment_urn,
        score_data=f"{test_file_dir}/test_score_data.csv",
        count_data=f"{test_file_dir}/test_count.csv",
        meta_data=f"{test_file_dir}/test_metadata.json",
        licence=Licence(short_name='CC BY 4.0'),
    
        sra_ids=['SRP109119'],
        pubmed_ids=['23035249'],
        doi_ids=['10.1038/s41467-019-11526-w'],
    )

Instantiate the NewScoresetRequest object and assign it to the
new_scoreset_request. You must substitute the attribute values for your
own.

.. code:: ipython3

    # substitute each attribute for your scoreset attributes
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

POST the model instance by passing the NewExperiment object as an
argument to the post_model_istance funtion that operates on the Client
object. This will POST the model instance to the approprate API
endpoint.

.. code:: ipython3

    client.post_model_instance(new_scoreset_request)

