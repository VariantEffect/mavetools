import asyncio
import logging
from .base import BaseClient
from aiohttp import ClientSession


class Client(BaseClient):
    """
    The Client object inherits the BaseClient and upon instantiation sets the base url where API requests will be made.
    CRUD operations can be made using the client object.
    """

    async def get_experiment(self, urn):
        """
        Hit an API endpoint to get instance of experiment by passing the experiment URN value.
        Parsed JSON data is returned.

        Parameters
        ----------
        urn : str
            The URN of the experiment to be retrieved.

        Returns
        -------
        str
            The experiment requested as a JSON string.
        """
        return await self.get_dataset("experiments", urn)

    async def get_experiments(self, urn_list):
        """
        Get a list of experiments by passing a list of URN values corresponding to the experiments.

        Parameters
        ----------
        urn_list: list

        Returns
        -------


        """
        async with ClientSession() as self.session:
            r = await asyncio.gather(*[self.get_dataset("experiments", urn) for urn in urn_list])
        return r

    async def get_scoreset(self, urn):
        """
        Hit an API endpoint to get instance of scoreset by passing the experiment URN value.
        Parsed JSON data is returned.

        Parameters
        ----------
        urn : str
            The URN of the scoreset to be retrieved.

        Returns
        -------
        str
            The scoreset requested as a JSON string.
        """
        async with ClientSession() as self.session:
            r = await self.get_dataset("scoresets", urn)
        return r

    async def get_scoresets(self, urn_list):
        """
        Get a list of scoresets by passing a list of URN values corresponding to the scoresets.

        Parameters
        ----------
        urn_list

        Returns
        -------

        """
        async with ClientSession() as self.session:
            r = await asyncio.gather(*[self.get_dataset("scoresets", urn) for urn in urn_list])
        return r

    async def create_experiment(self, experiment):
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
        return await self.create_dataset(experiment, "experiments")

    async def create_experiments(self, experiment_list):
        """
        Create more than one experiment by passing a list of experiments formatted as dictionaries.

        Parameters
        ----------
        experiment_list

        Returns
        -------

        """
        async with ClientSession() as self.session:
            r = await asyncio.gather(*[self.create_dataset(experiment, "experiments") for experiment in experiment_list])
        return r

    async def create_scoreset(self, scoreset, scores_df=None, counts_df=None):
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
        return await self.create_dataset(scoreset, "scoresets", scores_df, counts_df)

    async def create_scoresets(self, scoreset_list):
        """
        Create more than one scoreset by passing a list of scoresets formatted as dictionaries.

        Parameters
        ----------
        scoreset_list

        Returns
        -------

        """
        async with ClientSession() as self.session:
            r = await asyncio.gather(*[self.create_dataset(scoreset[0],
                                                           "scoresets",
                                                           scoreset[1],
                                                           scoreset[2]) for scoreset in scoreset_list])
        return r


