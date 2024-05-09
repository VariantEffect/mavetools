import os
import ssl
import certifi
import aiohttp
import pandas as pd
import humps
from urllib.parse import urlparse
from aiohttp import ClientResponseError
from mavedb.lib.validation.urn_re import (
    MAVEDB_SCORE_SET_URN_RE,
    MAVEDB_EXPERIMENT_URN_RE,
    MAVEDB_EXPERIMENT_SET_URN_RE,
)
from typing import Optional, Awaitable, Mapping
from mavetools.client.util import infer_record_type, validate_dataset_with_create_model
from mavedb.lib.validation.dataframe import validate_and_standardize_dataframe_pair
from mavedb.lib.validation.exceptions import ValidationError

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
            This token is required to enable data deposition (POST) operations and private record access.
        """
        if base_url is None:
            if os.environ.get(MAVEDB_API_URL) is None:
                raise ValueError(
                    f"API base URL not provided and not defined in OS environment under '{MAVEDB_API_URL}'"
                )
            else:
                base_url = os.environ.get(MAVEDB_API_URL)

        # split the base_url into the base and api_root portions since aiohttp doesn't allow paths in base_url
        parse_result = urlparse(base_url)
        if parse_result.scheme:
            base_url = f"{parse_result.scheme}://{parse_result.netloc}/"
        else:
            base_url = f"//{parse_result.netloc}/"
        self.api_root = parse_result.path

        self.session = aiohttp.ClientSession(
            base_url=base_url,
            connector=aiohttp.TCPConnector(ssl=ssl.create_default_context(cafile=certifi.where())),
            raise_for_status=True,
        )
        if auth_token is None:
            self.auth_token = ""
        else:
            self.auth_token = auth_token

        self.endpoints = {
            "score_set": "score-sets",
            "experiment": "experiments",
            "experiment_set": "experiment-sets",
        }

    async def get_dataset(self, urn: str, record_type: Optional[str] = None) -> Awaitable[str]:
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
            if MAVEDB_SCORE_SET_URN_RE.match(urn):
                record_type = "score_set"
            elif MAVEDB_EXPERIMENT_URN_RE.match(urn):
                record_type = "experiment"
            elif MAVEDB_EXPERIMENT_SET_URN_RE.match(urn):
                record_type = "experiment_set"
            else:
                raise ValueError(f"unable to infer record_type for '{urn}'")
        elif record_type not in self.endpoints.keys():
            raise ValueError(f"invalid record_type '{record_type}'")

        url_path = "/".join(x.strip("/") for x in ("", self.api_root, self.endpoints[record_type], urn))
        try:
            async with self.session.get(url_path, headers={"X-API-key": self.auth_token}) as resp:
                return await resp.json()
        except ClientResponseError as e:
            print(f"error {e.status} while requesting {url_path}")

    async def create_dataset(
        self, dataset: Mapping, scores_df: Optional[pd.DataFrame] = None, counts_df: Optional[pd.DataFrame] = None
    ) -> Optional[str]:
        """
        Submit a dataset to the API.

        If the dataset being submitted is a score set, it must include a `scores_df` and optional `counts_df`.

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
        Optional[str]
            The URN of the created dataset if successful; else None.

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
            r = await self.session.request(
                method="POST", url=url_path, json=dataset, headers={"X-API-key": self.auth_token}
            )
            r.raise_for_status()
            dataset = await r.json()
            dataset = humps.decamelize(dataset)
            urn = dataset["urn"]
        except ClientResponseError as e:
            print(f"error response {e.status} while requesting {url_path}")

        if record_type == "score_set" and urn is not None:
            await self.upload_dataframes(dataset, scores_df, counts_df)

        # return the URN of the created model instance
        return urn

    async def upload_dataframes(
        self, score_set: Mapping, scores_df: pd.DataFrame, counts_df: Optional[pd.DataFrame] = None
    ) -> None:
        """
        Validate and upload data frames for a score set.

        Parameters
        ----------
        score_set: Mapping
            Object for the score set associated with these data frames.
            Used to get the URN and target sequence information.
        scores_df: pd.DataFrame
            Pandas data frame containing the scores.
        counts_df: Optional[pd.DataFrame]
            Pandas data frame containing the counts (if available).

        Returns
        -------
        None

        Raises
        ------
        None

        """
        url_path = "/".join(
            x.strip("/") for x in ("", self.api_root, self.endpoints["score_set"], score_set["urn"], "variants", "data")
        )
        """
        # TODO: this needs to be updated for the current multi-target validator
        try:
            new_scores_df, new_counts_df = validate_and_standardize_dataframe_pair(
                scores_df,
                counts_df,
                score_set["target_genes"],
                None,
            )
        except ValidationError as e:
            print(f"data frames for '{score_set['urn']}' failed to validate: {e}")
            return
        """
        # TODO test this with a really big dataframe
        upload_data = dict()
        upload_data["scores_file"] = bytes(scores_df.to_csv(index=False), encoding="utf-8")
        if counts_df is not None:
            upload_data["counts_file"] = bytes(counts_df.to_csv(index=False), encoding="utf-8")

        try:  # to post data
            r = await self.session.request(
                method="POST", url=url_path, data=upload_data, headers={"X-API-key": self.auth_token}
            )
            r.raise_for_status()
        except ClientResponseError as e:
            print(f"error response {e.status} while uploading data to {url_path}")
