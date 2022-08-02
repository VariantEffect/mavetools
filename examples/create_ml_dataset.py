from mavetools.client.client import LocalClient
from mavetools.models.ml_tools import MlDataset

#This example shows how to create a scaled dataset of all SAV effect values contained in the MaveDB

#We use a locally cloned version here, if you haven't cloned it yet, look into clone.py

#Provide the path to the local clone
local_instance_path = f'../../localMaveDB'

#Provide paths to where the dataset will be written.
outfile = 'mave_db_scaled_savs.fasta'
seq_file = 'mave_db_scaled_savs_only_sequences.fasta'

#Create a local client object
client = LocalClient(local_instance_path)

#Search the database without any filters to retrieve the whole database
experiment_dict = client.search_database()

#Create the ML dataset object
ml_dataset = MlDataset(experiment_dict)

#Retrieve the scoretables
ml_dataset.retrieve_data(client)

#Aggregate scoresets from same experiments
ml_dataset.aggregate_scoresetdata()

#Scale all SAV effect scores
ml_dataset.scale_all_savs()

#Write the output
ml_dataset.write_scaled_sav_fasta(outfile, sequences_only_file = seq_file)
