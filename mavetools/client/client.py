from .base import BaseClient
import logging


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
        return self.get_dataset("experiments", urn)

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
        return self.get_dataset("scoresets", urn)

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
        # validate here
        return self.create_dataset(experiment, "experiments")

    def create_scoreset(self, scoreset, scores_df=None, counts_df=None):
        """
        Hit an API endpoint to post a scoreset.

        Parameters
        ----------
        scoreset: dict
            Instance of the scoreset that will be POSTed.
        scores_df: pandas.DataFrame
            The scores file associated with a ScoreSet.
        counts_df: pandas.DataFrame
            The counts file associated with a ScoreSet.

        Returns
        -------
        str
            The URN of the created model instance.

        Raises
        ------
        AuthTokenMissingException
            If the auth_token is missing
        ValueError
            If the dataset is a ScoreSet and there is no scores_df provided.
        """
        if scores_df is None:
            error_message = "Must include a scores_df when creating a ScoreSet!"
            logging.error(error_message)
            raise ValueError(error_message)
        # validate here
        return self.create_dataset(scoreset, "scoresets", scores_df, counts_df)




