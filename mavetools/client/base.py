import json
import logging
import requests
import httpx
import sys
import time


class BaseClient:
    """
    The BaseClient object upon instantiation sets the base url where API requests will be made.
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
        if auth_token:
            self.auth_token = auth_token

    class AuthTokenMissingException(Exception):
        pass

    def get_dataset(self, endpoint, urn):
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
            #time.sleep(0.01)
            #r = requests.get(instance_url)
            print("before request")
            print(instance_url)
            r = httpx.get(instance_url)
            print("after request")
            r.raise_for_status()
        #except requests.exceptions.HTTPError as e:
        except httpx.HTTPError as exc:
            print(f"Error response {exc.response.status_code} while requesting {exc.request.url!r}.")
            #print(logging.error(r.json()))
            #raise SystemExit(e)

        return r.json()

    def create_dataset(self, dataset, endpoint, scores_df=None, counts_df=None):
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
            print("Before request")
            print(model_url)
            r = requests.post(
                model_url,
                json=dataset,
                headers={"X-API-key": self.auth_token},
            )
            print("After request")
            r.raise_for_status()
            urn = json.loads(r.text)['urn']
        except requests.exceptions.HTTPError as e:
        #except httpx.HTTPError as exc:
            #print(f"Error response {exc.response.status_code} while requesting {exc.request.url!r}.")
            # see what this error looks like, present them in understandable way
            logging.error(r.text)
            #sys.exit(1)

        if scores_df is not None and urn is not None:
            model_url = f"{self.base_url}scoresets/{urn}/variants/data"
            file_upload = dict()
            file_upload["scores_file"] = bytes(scores_df.to_csv(), encoding='utf-8')
            if counts_df is not None: file_upload["counts_file"] = bytes(counts_df.to_csv(), encoding='utf-8')

            try:  # to post data
                print("Before request")
                print(model_url)
                r = httpx.post(
                    model_url,
                    files=file_upload,
                    headers={"X-API-key": self.auth_token},
                    timeout=None,
                )
                print("After request")
                r.raise_for_status()
            #except requests.exceptions.HTTPError as e:
            except httpx.HTTPError as exc:
                print(f"Error response {exc.response.status_code} while requesting {exc.request.url!r}.")
                #logging.error(r.text)
                #sys.exit(1)

        # No errors or exceptions at this point, log successful upload
        logging.info(f"Successfully uploaded {dataset}!")

        # return the URN of the created model instance
        return urn
