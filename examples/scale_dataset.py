from mavetools.client.client import LocalClient
from mavetools.models.ml_tools import MlDataset

#This example shows how to create a scaled for one particular scoreset

#We use a locally cloned version here, if you haven't cloned it yet, look into clone.py

#Give the path to local clone
local_instance_path = f'../../localMaveDB'

#Create the local client object
client = LocalClient(local_instance_path)

#Provide a MaveDB urn identifier
particular_experiment_id = 'urn:mavedb:00000005-a'

#Create ML dataset object and aggregate all scoresets
experiment_dict = client.get_experiment_dict([particular_experiment_id])

print('=== Experiment object created ===')

ml_dataset = MlDataset(experiment_dict)

print('=== ML object created ===')

ml_dataset.retrieve_data(client, verbosity = 1)

print('=== Data retrieval done ===')

ml_dataset.aggregate_scoresetdata(verbosity = 1)

print('=== Data aggregation done ===')

#Uncomment to plot the SAV score histogram
ml_dataset.experiments[particular_experiment_id].experiment_scoresetdata.plot_sav_score_distribution(f'score_distribution_{particular_experiment_id}.png')
ml_dataset.experiments[particular_experiment_id].experiment_scoresetdata.write_nonsense_tsv(f'Nonsense_scores_{particular_experiment_id}.tsv')

#Scale all SAV effect scores
ml_dataset.scale_all_savs(verbosity = 1)

print('=== Data scaling finished ===')

#Uncomment to plot the scaled SAV score histogram
ml_dataset.experiments[particular_experiment_id].experiment_scoresetdata.plot_sav_score_distribution(f'scaled_score_distribution_{particular_experiment_id}.png', specific_scores = ml_dataset.experiments[particular_experiment_id].experiment_scoresetdata.scaled_sav_scores)

#Write the scaled the dataset to the specialized fasta file
outfile = f'{particular_experiment_id}_scaled_savs.fasta'
ml_dataset.write_scaled_sav_fasta(outfile)
