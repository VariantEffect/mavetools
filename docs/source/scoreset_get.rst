GET scoreset
============

get a scoreset from MaveDB via the API
--------------------------------------

To begin, import the modeules below.

.. code:: ipython3

    import attr, os
    from pprint import PrettyPrinter
    from mavetools.client.client import Client
    from mavetools.models.scoreset import ScoreSet

Pretty printer is used to format the output nicely.

.. code:: ipython3

    pp = PrettyPrinter(indent=2)

Here your base_url is set to localhost, http://127.0.0.1:8000/api/. This
default funcionality is what you would want to use when working with a
local instance of MaveDB (e.g., a development branch). If working with
production mavedb you would set base url to https://www.mavedb.org/api/.

In the cell below, comment out the base_url you will not be using.

.. code:: ipython3

    base_url = 'http://127.0.0.1:8000/api/'
    #base_url = 'https://www.mavedb.org/api/'

Set experiment_urn to match the scoreset you want to get.

.. code:: ipython3

    scoreset_urn = 'urn:mavedb:00000001-a-1'

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

    # this is an example of what your auth_token should look like
    auth_token = 'R2skRbpBD3Rsf5dNHoQxDZevdEE74T5lCKMFyBhBwwPFH4ZfTrxDz7TZ0kbFLtEZ'

Here you instantiate the Client object. The Client object is the object
by which the POST request is performed. The client object is
instantiated with the value of base_url provided earlier, so make sure
that is up-to-date. If base_url does not exist, base_url is defaulted to
localhost, http://127.0.0.1:8000/api/.

.. code:: ipython3

    client = Client(base_url, auth_token=auth_token) if base_url else Client(auth_token=auth_token)

GET the model instance by passing the model type (Scoreset, in this
instance) and the scoreset_urn as arguments to the get_model_istance
funtion that operates on the Client object. This will GET the model
instance (resource) from the server via the approprate API endpoint.

.. code:: ipython3

    scoreset = client.get_model_instance(ScoreSet, scoreset_urn)

Now, display the results!

.. code:: ipython3

    pp.pprint(attr.asdict(scoreset))


.. parsed-literal::

    { 'abstract_text': 'Although we now routinely sequence human genomes, we can '
                       'confidently identify only a fraction of the sequence '
                       'variants that have a functional impact. Here, we developed '
                       'a deep mutational scanning framework that produces '
                       'exhaustive maps for human missense variants by combining '
                       'random codon mutagenesis and multiplexed functional '
                       'variation assays with computational imputation and '
                       'refinement. We applied this framework to four proteins '
                       'corresponding to six human genes: UBE2I (encoding SUMO E2 '
                       'conjugase), SUMO1 (small ubiquitin-like modifier), TPK1 '
                       '(thiamin pyrophosphokinase), and CALM1/2/3 (three genes '
                       'encoding the protein calmodulin). The resulting maps '
                       'recapitulate known protein features and confidently '
                       'identify pathogenic variation. Assays potentially amenable '
                       'to deep mutational scanning are already available for 57% '
                       'of human disease genes, suggesting that DMS could '
                       'ultimately map functional variation for all human disease '
                       'genes. \r\n'
                       '\r\n'
                       'See [**Weile *et al.* '
                       '2017**](http://msb.embopress.org/content/13/12/957)',
      'approved': None,
      'contributors': ['0000-0003-1628-9390'],
      'count_columns': ['hgvs_nt', 'hgvs_splice', 'hgvs_pro'],
      'created_by': '0000-0003-1628-9390',
      'creation_date': '2018-06-26',
      'current_version': 'urn:mavedb:00000001-a-1',
      'data_usage_policy': '',
      'dataset_columns': None,
      'doi_ids': [],
      'experiment': 'urn:mavedb:00000001-a',
      'extra_metadata': {},
      'is_meta_analysis': False,
      'keywords': [ {'text': 'DMS-BarSeq'},
                    {'text': 'E2'},
                    {'text': 'sumoylation'},
                    {'text': 'imputation'},
                    {'text': 'DMS-TileSeq'},
                    {'text': 'complementation'}],
      'last_child_value': None,
      'licence': { 'link': 'https://creativecommons.org/licenses/by/4.0/',
                   'long_name': 'CC BY 4.0 (Attribution)',
                   'short_name': 'CC BY 4.0',
                   'version': '4.0'},
      'method_text': '##Scoring procedure:\r\n'
                     'DMS-BarSeq and DMS-TileSeq reads were processed using the '
                     '[dmsPipeline](https://bitbucket.org/rothlabto/dmspipeline) '
                     'software. Briefly, Barseq read counts were used to establish '
                     'relative frequencies of each strain at each timepoint and '
                     'converted to estimates of absolute frequencies using OD '
                     'measurement data. Absolute counts were used to establish '
                     'growth curves from which fitness parameters were estimated '
                     'and then normalized to 0-1 scale where 0 corresponds to null '
                     'controls and 1 corresponds to WT controls. Meanwhile, '
                     'TileSeq read counts were used to establish relative allele '
                     'frequencies in each condition. Non-mutagenized control '
                     'counts were subtracted from counts (as estimates of '
                     'sequencing error). log ratios of selection over '
                     'non-selection counts were calculated. The resulting TileSeq '
                     'fitness values were then rescaled to the distribution of the '
                     'BarSeq fitness scores. Fitness scores were joined using '
                     'confidence-weighted averages. Random-Forest base machine '
                     'learning was used to impute missing values and refine '
                     'low-confidence measurements, based on intrinsic, structural, '
                     'and biochemical features.\r\n'
                     '\r\n'
                     'See [**Weile *et al.* '
                     '2017**](http://msb.embopress.org/content/13/12/957) for more '
                     'details.\r\n'
                     '\r\n'
                     '## Additional columns:\r\n'
                     '* exp.score = experimental score from the joint '
                     'DMS-BarSeq/DMS-TileSeq screens\r\n'
                     '* exp.sd = standard deviation of the experimental score\r\n'
                     '* df = degrees of freedom (number of replicates contributing '
                     'to the experimental score)\r\n'
                     '* pred.score = machine-learning predicted score',
      'modification_date': '2019-08-08',
      'modified_by': '0000-0003-1628-9390',
      'next_version': None,
      'previous_version': None,
      'private': None,
      'publish_date': '2018-06-26',
      'pubmed_ids': [ { 'dbname': 'PubMed',
                        'dbversion': None,
                        'identifier': '29269382',
                        'url': 'http://www.ncbi.nlm.nih.gov/pubmed/29269382'}],
      'replaces': None,
      'score_columns': [ 'hgvs_nt',
                         'hgvs_splice',
                         'hgvs_pro',
                         'score',
                         'sd',
                         'se',
                         'exp.score',
                         'exp.sd',
                         'df',
                         'pred.score'],
      'short_description': 'A joint Deep Mutational Scan of the human SUMO E2 '
                           'conjugase UBE2I using functional complementation in '
                           'yeast, combining DMS-BarSeq and DMS-TileSeq data, '
                           'followed by machine-learning-based imputation and '
                           'refinement.',
      'sra_ids': None,
      'target': { 'ensembl': { 'dbname': 'Ensembl',
                               'dbversion': None,
                               'identifier': 'ENSG00000103275',
                               'offset': 0,
                               'url': 'http://www.ensembl.org/id/ENSG00000103275'},
                  'name': 'UBE2I',
                  'reference_maps': [ { 'genome': { 'assembly_identifier': { 'dbname': 'GenomeAssembly',
                                                                             'dbversion': None,
                                                                             'identifier': 'GCF_000001405.26',
                                                                             'url': 'http://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.26'},
                                                    'organism_name': 'Homo sapiens',
                                                    'short_name': 'hg38'}}],
                  'reference_sequence': { 'sequence': 'ATGTCGGGGATCGCCCTCAGCAGACTCGCCCAGGAGAGGAAAGCATGGAGGAAAGACCACCCATTTGGTTTCGTGGCTGTCCCAACAAAAAATCCCGATGGCACGATGAACCTCATGAACTGGGAGTGCGCCATTCCAGGAAAGAAAGGGACTCCGTGGGAAGGAGGCTTGTTTAAACTACGGATGCTTTTCAAAGATGATTATCCATCTTCGCCACCAAAATGTAAATTCGAACCACCATTATTTCACCCGAATGTGTACCCTTCGGGGACAGTGTGCCTGTCCATCTTAGAGGAGGACAAGGACTGGAGGCCAGCCATCACAATCAAACAGATCCTATTAGGAATACAGGAACTTCTAAATGAACCAAATATCCAAGACCCAGCTCAAGCAGAGGCCTACACGATTTACTGCCAAAACAGAGTGGAGTACGAGAAAAGGGTCCGAGCACAAGCCAAGAAGTTTGCGCCCTCATAA',
                                          'sequence_type': 'dna'},
                  'refseq': { 'dbname': 'RefSeq',
                              'dbversion': None,
                              'identifier': 'NM_003345',
                              'offset': 159,
                              'url': 'http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?val=NM_003345'},
                  'scoreset': 'urn:mavedb:00000001-a-1',
                  'type': 'Protein coding',
                  'uniprot': { 'dbname': 'UniProt',
                               'dbversion': None,
                               'identifier': 'P63279',
                               'offset': 0,
                               'url': 'http://purl.uniprot.org/uniprot/P63279'}},
      'title': 'UBE2I imputed & refined',
      'urn': 'urn:mavedb:00000001-a-1',
      'variant_count': 3180}


