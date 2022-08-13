from .base import BaseClient


class Client(BaseClient):
    """
    The Client object inherits the BaseClient and upon instantiation sets the base url where API requests will be made.
    CRUD operations can be made using the client object.
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
        return self.get_model_instance("experiments", urn)

    def get_scoreset(self, urn):
        """
        Hit an API endpoint to get instance of scoreset by passing the experiment URN value.
        Parsed JSON data is returned.

        Parameters
        ----------
        urn : str
            The URN of the scoreset to be retrieved.

        Returns
        -------
        The scoreset requested.
        """
        return self.get_model_instance("scoresets", urn)

    def create_experiment(self, experiment):
        """
        Hit an API endpoint to post an experiment.

        Parameters
        ----------
        experiment: dict
            Instance of the experiment that will be POSTed.

        Returns
        -------
        str
            The URN of the created model instance.

        Raises
        ------
        AuthTokenMissingException
            If the auth_token is missing
        """
        return self.post_model_instance(experiment)

    def create_scoreset(self, scoreset):
        """
        Hit an API endpoint to post a scoreset.

        Parameters
        ----------
        scoreset: dict
            Instance of the scoreset that will be POSTed.
        files:
            The files associated with the scoreset.

        Returns
        -------
        str
            The URN of the created model instance.

        Raises
        ------
        AuthTokenMissingException
            If the auth_token is missing
        """
        return self.post_model_instance(scoreset)




