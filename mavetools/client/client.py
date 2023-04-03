import asyncio
import logging
import ssl
import certifi
from aiohttp import ClientSession, TCPConnector, ClientResponseError


class Client:
    """
    The Client object sets the base url where API requests will be made.
    CRUD operations can be made using the client object.
    """
    def __init__(self, base_url="http://127.0.0.1:8002/api/v1/", auth_token=""):
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
        self.session = None
        self.sslcontext = ssl.create_default_context(cafile=certifi.where())
        if auth_token:
            self.auth_token = auth_token

    class AuthTokenMissingException(Exception):
        pass

    async def get_dataset(self, endpoint, urn):
        """
        Using a GET, hit an API endpoint to get a dataset such as a ScoreSet.
        This will perform the HTTP GET request and return the dataset as a
        JSON string.

        Parameters
        ----------
        endpoint : str
            The API endpoint where we want the request to be made. This is the url extension beyond the base url
            used to instantiate the Client object. For example if you want an experiment from the base url
            'http://127.0.0.1:8000/api/v1/', the api_endpoint argument would be "experiments", making an API endpoint
            of 'http://127.0.0.1:8000/api/v1/experiments'.
        urn : str
            The URN of the object we are retrieving.

        Returns
        -------
        str
            An instance of the dataset as a JSON str.

        Raises
        ------
        ValueError
            If any mandatory fields are missing.
        """
        model_url = f"{self.base_url}{endpoint}/"
        instance_url = f"{model_url}{urn}"
        try:
            r = await self.session.request(method="GET", url=instance_url)
            r.raise_for_status()
        except ClientResponseError as e:
            print(f"Error response {e.status} while requesting {instance_url!r}.")

        return await r.json()

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
        async with ClientSession(connector=TCPConnector(ssl=self.sslcontext)) as self.session:
            r = await asyncio.gather(*[self.create_dataset(scoreset[0],
                                                           "scoresets",
                                                           scoreset[1],
                                                           scoreset[2]) for scoreset in scoreset_list])
        return r


