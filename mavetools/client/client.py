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

    def get_experiment(self, urn):
        """
        Hit an API endpoint to get instance of experiment by passing the experiment URN value.
        Parsed JSON data is returned.

        Parameters
        ----------
        urn : str
            The URN of the experiment to be retrieved.

        Returns
        -------
        The experiment requested.
        """
        return self.get_model_instance(Experiment, urn)

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




