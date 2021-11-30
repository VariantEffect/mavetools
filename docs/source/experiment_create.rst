POST experiment
===============

create an experiment in MaveDB via the API
------------------------------------------

To begin, import the modeules below.

.. code:: ipython3

    import os
    from mavetools.client.client import Client
    from mavetools.models.experiment import NewExperiment

Here your base_url is set to localhost, http://127.0.0.1:8000/api/. This
default funcionality is what you would want to use when working with a
local instance of MaveDB (e.g., a development branch). If working with
production mavedb you would set base url to https://www.mavedb.org/api/.

In the cell below, comment out the base_url you will not be using.

.. code:: ipython3

    base_url = 'http://127.0.0.1:8000/api/'
    #base_url = 'https://www.mavedb.org/api/'

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

For the abstrct.md and method.md files, simply paste your content into
these files.

We have an example directory in mavetools that holds the files of
interest. Though this directory can exist anywhere on your computer, you
must put the correct path to that directory as the value to
test_file_dir.

.. code:: ipython3

    # here is an example if your copy of mavetools exists within Pycharm
    test_file_dir = '~/PycharmProjects/mavetools/tests/test_upload_scoreset/test_files'

Here, we will want to validate the data that will be POSTed. Run the
validation code.

.. code:: ipython3

    #validate

Create a NewExperiment object and populate the required attributes. You
will need to substitute the values of the attributes for your own.

.. code:: ipython3

    #with open(f"{test_file_dir}/abstract.md") as handle:
    #    test_abstract_text = handle.read()
    #with open(f"{test_file_dir}/method.md") as handle:
    #    test_method_text = handle.read()
        
    # substitute the each attribute for your experiment attributes
    new_experiment = NewExperiment(
        # experimentset=experimentset_urn,
        title='exp_test_title',
        short_description='exp_test_short_description',
        #abstract_text=test_abstract_text,
        #method_text=test_method_text,
        abstract_text="abstract",
        method_text="method",
    
        sra_ids=['SRP109119'],
        pubmed_ids=['23035249'],
        doi_ids=['10.1038/s41467-019-11526-w'],
    )

POST the model instance by passing the NewExperiment object as an
argument to the post_model_istance funtion that operates on the Client
object. This will POST the model instance to the approprate API
endpoint.

.. code:: ipython3

    client.post_model_instance(new_experiment)

