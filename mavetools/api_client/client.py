import json
import requests
import logging
import sys


class Client:
    def __init__(self, base_url="http://127.0.0.1:8000/api/", auth_token=""):
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

    def get_experiment(self, instance_id):
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
        model_url = f"{self.base_url}v1/experiments/"
        instance_url = f"{model_url}{instance_id}/"
        try:
            r = requests.get(instance_url)
            r.raise_for_status()
        except requests.exceptions.HTTPError as e:
            logging.error(r.json())
            raise SystemExit(e)
        return r.json()

    def get_scoreset(self, instance_id):
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
        model_url = f"{self.base_url}v1/scoresets/"
        instance_url = f"{model_url}{instance_id}/"
        try:
            r = requests.get(instance_url)
            r.raise_for_status()
        except requests.exceptions.HTTPError as e:
            logging.error(r.json())
            raise SystemExit(e)
        print(r)
        print(r.json())
        return r.json()

    def post_experiment(self, model_instance):
        """
        Using a POST, hit an API endpoint to post a resource.
        Performs HTTP POST request.

        Parameters
        ----------
        model_instance
            instance of model that will be POSTed

        Returns
        -------
        AuthTokenMissingException
            If the auth_token is missing
        """

        # save object type of model_instance
        #model_class = type(model_instance)
        #model_url = f"{self.base_url}{model_class.api_url()}/"
        model_url = f"{self.base_url}v1/experiments/"
        payload, files = model_instance.post_payload()

        # check for existence of self.auth_token, raise error if does not exist
        if not self.auth_token:
            error_message = "Need to include an auth token for POST requests!"
            logging.error(error_message)
            raise self.AuthTokenMissingException(error_message)

        try:  # to post data
            r = requests.post(
                model_url,
                #data=model_instance,
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