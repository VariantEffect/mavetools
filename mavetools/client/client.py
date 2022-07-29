import json
import logging
import requests
import sys
import os

from mavetools.models.scoreset import ScoreSet
from mavetools.models.ml_tools import MlExperiment

class ClientTemplate:
    def parse_json_scoreset_list(self, scoreset_list, keywords = None, organisms = None, retrieve_json_only = False, experiment_types = None):
        if keywords is not None:
            keywords = set(keywords)

        if organisms is not None:
            organisms = set(organisms)

        if experiment_types is not None:
            experiment_types = set(experiment_types)

        experiment_dict = {}
        for scoreset in scoreset_list:
            if keywords is not None:
                keyword_match = False
                for keyword in scoreset['keywords']:
                    if keyword['text'] in keywords:
                        keyword_match = True
                if not keyword_match:
                    continue

            if organisms is not None:
                if not scoreset['target']['reference_maps'][0]['genome']['organism_name'] in organisms:
                    continue

            if experiment_types is not None:
                if not scoreset['target']['type'] in experiment_types:
                    continue

            urn = scoreset['urn']

            if retrieve_json_only:
                experiment_dict[urn] = scoreset
                continue

            scoreset_obj = ScoreSet.deserialize(scoreset)

            experiment_urn = scoreset_obj.experiment

            if not experiment_urn in experiment_dict:
                experiment_dict[experiment_urn] = MlExperiment(experiment_urn, {}, scoreset_obj)

            experiment_dict[experiment_urn].scoreset_dict[urn] = scoreset_obj

        return experiment_dict

class LocalClient(ClientTemplate):
    def __init__(self, local_instance_path):
        self.local_instance_path = local_instance_path
        self.meta_data_folder = f'{local_instance_path}/meta_data/'
        self.scoreset_data_folder = f'{local_instance_path}/scoreset_data/'

    def get_meta_file_path(self, urn):
        return f'{self.meta_data_folder}/{urn}.json'

    def load_meta_data(self, filepath):
        f = open(filepath, 'r')
        meta_data = json.load(f)
        f.close()
        return meta_data

    def get_meta_data(self, urn):
        return self.load_meta_data(self.get_meta_file_path(urn))

    def search_database(self, keywords = None, organisms = None, experiment_types = ['Protein coding']):
        scoreset_list = []
        for meta_data_file in os.listdir(self.meta_data_folder):
            meta_data_file = f'{self.meta_data_folder}{meta_data_file}'

            meta_data = self.load_meta_data(meta_data_file)

            scoreset_list.append(meta_data)

        experiment_dict = self.parse_json_scoreset_list(scoreset_list, keywords = keywords, organisms = organisms, experiment_types = experiment_types)

        return experiment_dict

    def get_experiment_dict(self, urns):
        scoreset_list = []
        for urn in urns:
            scoreset_list.append(self.get_meta_data(urn))
        return self.parse_json_scoreset_list(scoreset_list)

    def retrieve_score_table(self, urn):
        score_table_file = f'{self.scoreset_data_folder}/{urn}.csv'
        f = open(score_table_file, 'r')
        text = f.read()
        f.close()
        return text


class Client(ClientTemplate):
    def __init__(self, base_url="https://www.mavedb.org/api/", auth_token=""):
        """
        Instantiates the Client object and sets the values for base_url and
        auth_token

        Parameters
        ----------
        base_url: the url in which the api endpoint exists
            default: 'http://127.0.0.1:8000/api/'
        auth_token: authorizes POST requests via the API and MaveDB
            default: ''
        """
        self.base_url = base_url
        if auth_token:
            self.auth_token = auth_token

    class AuthTokenMissingException(Exception):
        pass

    def clone(self, local_instance_path):
        if not os.path.exists(local_instance_path):
            os.mkdir(local_instance_path)
        meta_data_folder = f'{local_instance_path}/meta_data/'
        if not os.path.exists(meta_data_folder):
            os.mkdir(meta_data_folder)
        scoreset_data_folder = f'{local_instance_path}/scoreset_data/'
        if not os.path.exists(scoreset_data_folder):
            os.mkdir(scoreset_data_folder)

        entry_dict = self.search_database(retrieve_json_only = True)
        for urn in entry_dict:
            meta_file = f'{meta_data_folder}/{urn}.json'

            f = open(meta_file,'w')
            json.dump(entry_dict[urn], f)
            f.close()

            score_table_file = f'{scoreset_data_folder}/{urn}.csv'

            f = open(score_table_file, 'w')
            f.write(self.retrieve_score_table(urn))
            f.close()

    def search_database(self, keywords = None, organisms = None, retrieve_json_only = False, experiment_types = ['Protein coding']):
        search_page_url = f'{self.base_url}/scoresets'
        try:
            r = requests.get(search_page_url)
            r.raise_for_status()
        except requests.exceptions.HTTPError as e:
            logging.error(r.json())
            raise SystemExit(e)

        scoreset_list = r.json()

        experiment_dict = self.parse_json_scoreset_list(scoreset_list, keywords = keywords, organisms = organisms, retrieve_json_only = retrieve_json_only)

        return experiment_dict

    def retrieve_score_table(self, urn):
        base_parent = self.base_url.replace('api/','')
        score_table_url = f'{base_parent}scoreset/{urn}/scores/'
        try:
            r = requests.get(score_table_url)
            r.raise_for_status()
        except requests.exceptions.HTTPError as e:
            logging.error(r.json())
            raise SystemExit(e)
        return r.text


    def get_model_instance(self, model_class, instance_id):
        """
        Using a GET, hit an API endpoint to get info on a particular instance
        of a model class such as a ScoreSet.
        This will perform the HTTP GET request and then let the class itself
        parse the JSON data.

        Parameters
        ----------
        model_class : ModelClass
            The model class we want to which we want to cast the response.
            (e.g., Experiment or Scoreset)
        instance_id : str
            The id of the object we are retrieving.

        Returns
        -------
        model_instance
            An instance of the passed class.

        Raises
        ------
        ValueError
            If any mandatory fields are missing.
        """
        model_url = f"{self.base_url}{model_class.api_url()}"
        instance_url = f"{model_url}{instance_id}/"
        try:
            r = requests.get(instance_url)
            r.raise_for_status()
        except requests.exceptions.HTTPError as e:
            logging.error(r.json())
            raise SystemExit(e)
        return model_class.deserialize(r.json())

    def post_model_instance(self, model_instance):
        """
        Using a POST, hit an API endpoint to post a resource.
        Performs HTTP POST request.

        Parameters
        ----------
        model_instance
            instance of model that will be POSTed

        Returns
        -------
        requests.model.Response
            The HTTP response object from the request, which contains the URN
            of the newly-created model in the `Response.text` field.

        Raises
        ------
        AuthTokenMissingException
            If the auth_token is missing
        """

        # save object type of model_instance
        model_class = type(model_instance)
        model_url = f"{self.base_url}{model_class.api_url()}/"
        payload, files = model_instance.post_payload()

        # check for existance of self.auth_token, raise error if does not exist
        if not self.auth_token:
            error_message = "Need to include an auth token for POST requests!"
            logging.error(error_message)
            raise AuthTokenMissingException(error_message)

        try:  # to post data
            r = requests.post(
                model_url,
                data={"request": json.dumps(payload)},
                files=files,
                headers={"Authorization": (self.auth_token)},
            )
            r.raise_for_status()
        except requests.exceptions.HTTPError as e:
            logging.error(r.text)
            sys.exit(1)

        # No errors or exceptions at this point, log successful upload
        logging.info(f"Successfully uploaded {model_instance}!")

        # return the HTTP response
        return r
