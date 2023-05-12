import os
import asyncio
import logging
import ssl
import certifi
import aiohttp
from urllib.parse import urlparse
from aiohttp import ClientResponseError
from mavedb.lib.validation.constants.urn import MAVEDB_SCORESET_URN_RE, MAVEDB_EXPERIMENT_URN_RE, MAVEDB_EXPERIMENTSET_URN_RE
from typing import Optional, Awaitable, Mapping
import pandas as pd
from mavetools.client.util import infer_record_type, validate_dataset_with_create_model

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

        self.endpoints = {
            "score_set" : "scoresets",
            "experiment" : "experiments",
            "experiment_set" : "experimentSets",
        }

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
        ValueError
            If the record_type is invalid.
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
        elif record_type not in self.endpoints.keys():
            raise ValueError(f"invalid record_type '{record_type}'")

        url_path = "/".join(x.strip("/") for x in ("", self.api_root, self.endpoints[record_type], urn))
        try:
            async with self.session.get(url_path) as resp:
                return await resp.json()
        except ClientResponseError as e:
            print(f"error {e.status} while requesting {url_path}")

    async def create_dataset(self, dataset : Mapping, scores_df : Optional[pd.DataFrame] = None, counts_df : Optional[pd.DataFrame] = None):
        """
        Submit a dataset to the API.

        If the dataset being submitted is a score set, it must include a `scores_df` and optional `counts_df`.

        #TODO: support other types for the scores and counts, such as .csv files already on disk

        Parameters
        ----------
        dataset: Mapping
            Instance of the dataset that will be POSTed.
        scores_df: Optional[pd.DataFrame]
            The scores file associated with score set.
        counts_df: Optional[pd.DataFrame]
            The counts file associated with score set.

        Returns
        -------
        str
            The URN of the created dataset.

        Raises
        ------
        ValueError
            If the auth_token is missing.
        ValueError
            If the dataset is a ScoreSet and there is no scores_df provided.
        """
        # check for existence of self.auth_token, raise error if does not exist
        if not self.auth_token:
            error_message = "client must have an auth token to create datasets"
            raise ValueError(error_message)

        record_type = infer_record_type(dataset)
        if record_type is None:
            raise ValueError("could not infer record type for dataset")

        if record_type == "score_set" and scores_df is None:
            raise ValueError("must include a scores_df when creating a score set")

        # perform validation
        validate_dataset_with_create_model(dataset)

        url_path = "/".join(x.strip("/") for x in ("", self.api_root, self.endpoints[record_type]))
        urn = None

        try:  # to post data
            r = await self.session.request(method="POST",
                                           url=url_path,
                                           json=dataset,
                                           headers={"X-API-key": self.auth_token})
            r.raise_for_status()
            dataset = await r.json()
            urn = dataset['urn']
        except ClientResponseError as e:
            print(f"error response {e.status} while requesting {url_path}")

        if record_type == "score_set":
            url_path = "/".join(x.strip("/") for x in ("", self.api_root, self.endpoints[record_type], urn, "variants", "data"))

            # TODO test this with a really big dataframe
            upload_data = dict()
            upload_data["scores_file"] = bytes(scores_df.to_csv(), encoding='utf-8')
            if counts_df is not None:
                upload_data["counts_file"] = bytes(counts_df.to_csv(), encoding='utf-8')

            try:  # to post data
                r = await self.session.request(method="POST",
                                               url=url_path,
                                               data=upload_data,
                                               headers={"X-API-key": self.auth_token})
                r.raise_for_status()
            except ClientResponseError as e:
                print(f"error response {e.status} while uploading data to {url_path}")

        # return the URN of the created model instance
        return urn
