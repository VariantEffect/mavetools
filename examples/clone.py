import sys

from mavetools.client.client import Client

#This example shows how to download the whole MaveDB and create a local clone

#Provide the URL of MaveDB
base_url = 'https://www.mavedb.org/api/'

# Generate a new auth_token in your profile and post it here

auth_token = ""
# if the base url exists, the client object is instantiated with that value
# otherwise the client object is instantiated with default value which points to localhost
client = (
    Client(base_url, auth_token=auth_token)
    if base_url
    else Client(auth_token=auth_token)
)

#Provide a path where the local clone should be stored.
local_instance_path = f'../../localMaveDB'

#Download MaveDB
experiment_dict = client.clone(local_instance_path)

