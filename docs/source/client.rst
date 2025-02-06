Client
======
The Client module provides a Python interface for interacting with the MaveDB API.

Example usage
-------------

This page includes some examples of using the Client.

If you are looking to download a very large amount of data from MaveDB,
we strongly recommend using the versioned archive `available on Zenodo <https://doi.org/10.5281/zenodo.11201736>`_.

First, we'll set up a MaveDB Client instance using our API key, 
which can be generated on the MaveDB profile page::

   from mavetools.client.client import Client

   if "MAVEDB_API_KEY" in os.environ:
      api_key = os.environ.get("MAVEDB_API_KEY")
   else:
      api_key = "aaa-xxx-000"

   # API URL for local MaveDB instance
   # api_url = "http://localhost:8002/api/v1/"

   # API URL for MaveDB staging/testing instance
   # Note this may be running different version than the production server
   # api_url = "http://api.staging.mavedb.org/api/v1/"

   # API URL for the production MaveDB instance
   api_url = "https://api.mavedb.org/api/v1/"

   my_client = Client(base_url=api_url, auth_token=api_key)

We can then use the Client to retrieve a dataset::
      
   import pprint

   my_data = await my_client.get_dataset("urn:mavedb:00000013-a")
   pprint.pp(my_data)

The API can also be used to deposit datasets::

   import pandas as pd

   my_experiment = {
      "title": "Great Dataset",
      "short_description": "Very cool dataset where I did a neat experiment.",
      "abstract_text": "This is the abstract for my extremely cool experiment. It should have 2-3 sentences describing the study."
      "method_text": "Here are a few sentences summarizing the experimental methods I used. Since this is an experiment record, it does not include data analysis."
      "extra_metadata": {},
      "primary_publication_identifiers": [],
      "raw_read_identifiers": [],
   }
   new_experiment_urn = await my_client.create_dataset(my_experiment)
   print(f"deposited new experiment {new_experiment_urn}")

   my_score_set = {
      "experiment_urn": new_experiment_urn,
      "title": "Great Dataset Scores",
      "short_description": "Scores for the very cool dataset where I did a neat experiment.",
      "abstract_text": "This is the abstract for my extremely cool experiment. It should have 2-3 sentences describing the study. It is often the same as the experiment."
      "method_text": "Here are a few sentences summarizing the analysis methods I used. Since this is a score set record, it should start from the FASTQ files generated."
      "extra_metadata": {},
      "primary_publication_identifiers": [],
      "license_id": 1,
      "target_genes": [
         "name": "Target Gene",
         "category": "Protein coding",
         "external_identifiers": [
            {
               "identifier": {
                  "dbName": "UniProt",
                  "identifier": "UABC1234"
               },
               "offset": 12,
         ],
         "target_sequence": {
            "sequence": "ACGTTTACGTGG",
            "sequence_type": "dna",
            "taxonomy": {
               "tax_id": 9606,
            }
         }
      ],
   }
   new_score_set_urn = await my_client.create_dataset(my_score_set, 
                                                      scores_df=pd.read_csv(f"mavedb_files/scores.csv"),
                                                      counts_df=pd.read_csv(f"mavedb_files/counts.csv"),
                                                      )
   print(f"deposited new score set {new_score_set_urn}")

Note that we need to create an experiment record first, then include the MaveDB urn of that experiment in the score set.

After uploading a dataset, it will be stored with a temporary status and be visible only to you.
You should always log into the MaveDB web interface and inspect it before making it public.

API reference
-------------

.. automodule:: mavetools.client.client
   :members:

.. automodule:: mavetools.client.util
   :members:
