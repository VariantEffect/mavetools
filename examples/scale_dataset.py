from mavetools.client.client import LocalClient
from mavetools.models.ml_tools import MlDataset

#This example shows how to create a scaled for one particular scoreset

#We use a locally cloned version here, if you haven't cloned it yet, look into clone.py

#Give the path to local clone
local_instance_path = f'../../localMaveDB'

#Create the local client object
client = LocalClient(local_instance_path)

#Provide a MaveDB urn identifier
particular_experiment_id = 'urn:mavedb:00000045-c-1'

#Create ML dataset object and aggregate all scoresets
experiment_dict = client.get_experiment_dict([particular_experiment_id])
ml_dataset = MlDataset(experiment_dict)
ml_dataset.retrieve_data(client)
ml_dataset.aggregate_scoresetdata()

#Uncomment to plot the SAV score histogram
#ml_dataset.experiments[particular_experiment_id].experiment_scoresetdata.plot_sav_score_distribution('score_distribution.png')

#Scale all SAV effect scores
ml_dataset.scale_all_savs(verbosity = 1)

#Write the scaled the dataset to the specialized fasta file
outfile = f'{particular_experiment_id}_scaled_savs.fasta'
ml_dataset.write_scaled_sav_fasta(outfile)
