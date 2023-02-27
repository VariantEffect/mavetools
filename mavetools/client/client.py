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
        requests.model.Response
            The HTTP response object from the request, which contains the URN
            of the newly-created model in the `Response.text` field.

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
