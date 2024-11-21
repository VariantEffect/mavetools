import sys

from mavetools.client.client import LocalClient
from mavetools.models.ml_tools import MlDataset

if len(sys.argv < 3):
    print('Usage: python make_gold_standard.py [path to local instance of MaveDB (https://zenodo.org/records/11201737)] [path to the ProteinGym Substitutions reference file (https://marks.hms.harvard.edu/proteingym/DMS_ProteinGym_substitutions.zip)]')
    sys.exit(1)
#This example shows how to create a scaled dataset of all SAV effect values contained in the MaveDB and ProteiGym

#Provide the path to the local clone (download via https://zenodo.org/records/11201737)
local_instance_path = sys.argv[1]

#Provide the path to the ProteinGym Substitutions reference file (download via https://marks.hms.harvard.edu/proteingym/DMS_ProteinGym_substitutions.zip)
path_to_reference_file = sys.argv[2]

#Provide paths to where the dataset will be written.
outfile = 'mave_db_gold_standard.fasta'
seq_file = 'mave_db_gold_standard_only_sequences.fasta'

#Create a local client object
client = LocalClient(local_instance_path)

#Search the database without any filters to retrieve the whole database
experiment_dict = client.search_database()

#Create the ML dataset object
ml_dataset = MlDataset(experiment_dict)

#Retrieve the scoretables
ml_dataset.retrieve_data(client, verbosity = 1)

#Add the ProteinGym database
ml_dataset.load_protein_gym(path_to_reference_file)

#Aggregate scoresets from same experiments
filtered_list = ml_dataset.aggregate_scoresetdata(min_prot_size = 50, min_len_coverage = 0.4, std_filter = 0.25, verbosity = 1)

#ml_dataset.write_filtered_entries(filtered_list, 'filtered_entries_mave_db_gold_standard.tsv')

#Scale all SAV effect scores
ml_dataset.scale_all_savs(verbosity = 1)

#ml_dataset.plot_sav_score_distribution('mave_db_gold_standard_score_distribution.png')

#Write the output
ml_dataset.write_scaled_sav_fasta(outfile, sequences_only_file = seq_file)
