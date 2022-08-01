import json
import logging
import requests
import sys

from .base import BaseClient
from mavetools.models.experiment import Experiment
from mavetools.models.scoreset import ScoreSet


class Client(BaseClient):
    """

    """

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
        return self.get_model_instance(Experiment, instance_id)

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
        return self.get_model_instance(ScoreSet, instance_id)


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
        requests.model.Response
            The HTTP response object from the request, which contains the URN
            of the newly-created model in the `Response.text` field.

        Raises
        ------
        AuthTokenMissingException
            If the auth_token is missing
        """
        return self.post_model_instance(model_instance)

    def post_scoreset(self, model_instance):
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
        return self.post_model_instance(model_instance)

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
