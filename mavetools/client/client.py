import os
import asyncio
import logging
import ssl
import certifi
import aiohttp
from urllib.parse import urlparse
from aiohttp import ClientResponseError
from mavedb.view_models.scoreset import ScoresetCreate
from mavedb.view_models.experiment import ExperimentCreate
from mavedb.lib.validation.constants.urn import MAVEDB_SCORESET_URN_RE, MAVEDB_EXPERIMENT_URN_RE, MAVEDB_EXPERIMENTSET_URN_RE
from typing import Optional, Awaitable

MAVEDB_API_URL = "MAVEDB_API_URL"


class Client:
    """
    Client objects provide an object-oriented Python interface for the MaveDB API.
    """
    def __init__(self, base_url: Optional[str] = None, auth_token: Optional[str] = None):
        """
        Instantiate a new Client object.

        Parameters
        ----------
        base_url : Optional[str]
            The url for the server record_type, e.g. 'http://localhost:8002/api/v1/'.
            If this is None, the program will try to load an environment variable defined under MAVEDB_API_URL,
            defined in this file.
            If the url contains any queries or fragments, those will be discarded.
        auth_token: Optional[str]
            The API authorization token from the user's profile on the MaveDB server.
            This token is required to enable data deposition (POST) operations.
        """
        if base_url is None:
            if os.environ.get(MAVEDB_API_URL) is None:
                raise ValueError(f"API base URL not provided and not defined in OS environment under '{MAVEDB_API_URL}'")
            else:
                base_url = os.environ.get(MAVEDB_API_URL)

        # split the base_url into the base and api_root portions since aiohttp doesn't allow paths in base_url
        parse_result = urlparse(base_url)
        if parse_result.scheme:
            base_url = f"{parse_result.scheme}://{parse_result.netloc}/"
        else:
            base_url = f"//{parse_result.netloc}/"
        self.api_root = parse_result.path

        self.session = aiohttp.ClientSession(base_url=base_url, connector=aiohttp.TCPConnector(ssl=ssl.create_default_context(cafile=certifi.where())), raise_for_status=True)
        self.auth_token = auth_token

    class AuthTokenMissingException(Exception):
        pass

    async def get_dataset(self, urn : str, record_type : Optional[str] = None) -> Awaitable[str]:
        """
        Request a dataset from the API in JSON format.

        Parameters
        ----------
        urn : str
            The URN of the dataset being requested.
        record_type : Optional[str]
            The type of record to get, one of "score_set", "experiment", or "experiment_set".
            If this is None, the record_type will be inferred from the URN.
            Note that this must be provided for `tmp:` records.

        Returns
        -------
        Awaitable[str]
            The dataset in JSON format.

        Raises
        ------
        ValueError
            If the URN cannot be inferred.
        """
        # infer record_type if needed
        if record_type is None:
            if MAVEDB_SCORESET_URN_RE.match(urn):
                record_type = "score_set"
            elif MAVEDB_EXPERIMENT_URN_RE.match(urn):
                record_type = "experiment"
            elif MAVEDB_EXPERIMENTSET_URN_RE.match(urn):
                record_type = "experiment_set"
            else:
                raise ValueError(f"unable to infer record_type for '{urn}'")

        # set API endpoint
        if record_type == "score_set":
            endpoint = "scoresets"
        elif record_type == "experiment":
            endpoint = "experiments"
        elif record_type == "experiment_set":
            endpoint = "experimentSets"
        else:
            raise ValueError(f"invalid record_type '{record_type}'")

        url_path = "/".join(x.strip("/") for x in ("", self.api_root, endpoint, urn))
        try:
            async with self.session.get(url_path) as resp:
                return await resp.json()
        except ClientResponseError as e:
            print(f"error {e.status} while requesting {url_path}")

    async def create_dataset(self, dataset, endpoint, scores_df=None, counts_df=None):
        """
        Using an HTTP POST request, hit an API endpoint to create a dataset. When creating a Scoreset,
        you must include a scores_df to have a complete upload.

        Parameters
        ----------
        dataset: dict
            Instance of the dataset that will be POSTed.
        endpoint: str
            The API endpoint where we want the request to be made. This is the url extension beyond the base url
            used to instantiate the Client object. For example if you want an experiment from the base url
            'http://127.0.0.1:8000/api/v1/', the api_endpoint argument would be "experiments", making an API endpoint
            of 'http://127.0.0.1:8000/api/v1/experiments'.
        scores_df: pandas.DataFrame
            The scores file associated with a ScoreSet.
        counts_df: pandas.DataFrame
            The counts file associated with a ScoreSet.

        Returns
        -------
        str
            The URN of the created dataset.

        Raises
        ------
        AuthTokenMissingException
            If the auth_token is missing.
        ValueError
            If the dataset is a ScoreSet and there is no scores_df provided.
        """
        model_url = f"{self.base_url}{endpoint}/"
        urn = None

        # check for existence of self.auth_token, raise error if does not exist
        if not self.auth_token:
            error_message = "Must include an auth token when creating datasets!"
            logging.error(error_message)
            raise self.AuthTokenMissingException(error_message)

        if scores_df is None and endpoint == "scoresets":
            error_message = "Must include a scores_df when creating a ScoreSet!"
            logging.error(error_message)
            raise ValueError(error_message)

        try:  # to post data
            r = await self.session.request(method="POST",
                                           url=model_url,
                                           json=dataset,
                                           headers={"X-API-key": self.auth_token})
            r.raise_for_status()
            dataset = await r.json()
            urn = dataset['urn']
        except ClientResponseError as e:
            print(f"Error response {e.status} while requesting {model_url!r}.")

        if scores_df is not None and urn is not None:
            model_url = f"{self.base_url}scoresets/{urn}/variants/data"
            file_upload = dict()
            # TODO test this with a really big dataframe
            file_upload["scores_file"] = bytes(scores_df.to_csv(), encoding='utf-8')
            if counts_df is not None: file_upload["counts_file"] = bytes(counts_df.to_csv(), encoding='utf-8')

            try:  # to post data
                r = await self.session.request(method="POST",
                                               url=model_url,
                                               data=file_upload,
                                               headers={"X-API-key": self.auth_token})
                r.raise_for_status()
            except ClientResponseError as e:
                print(f"Error response {e.status} while requesting {model_url!r}.")

        # No errors or exceptions at this point, log successful upload
        logging.info(f"Successfully uploaded {dataset}!")

        # return the URN of the created model instance
        return urn

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
        async with ClientSession(connector=TCPConnector(ssl=self.sslcontext)) as self.session:
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
        async with ClientSession(connector=TCPConnector(ssl=self.sslcontext)) as self.session:  #connector=TCPConnector(ssl=False)
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
        async with ClientSession(connector=TCPConnector(ssl=self.sslcontext)) as self.session:
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
            error_message = "cannot create a new score set without a scores dataframe"
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
        async with ClientSession(connector=TCPConnector(ssl=self.sslcontext)) as self.session:
            r = await asyncio.gather(*[self.create_dataset(scoreset[0],
                                                           "scoresets",
                                                           scoreset[1],
                                                           scoreset[2]) for scoreset in scoreset_list])
        return r


