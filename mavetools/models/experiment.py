import attr
from typing import List

from .base import APIObject
from .dataset import Dataset, NewDataset
from .utils import attrs_filter, attrs_serializer, prepare_for_encoding


@attr.s
class Experiment(APIObject, Dataset):
    experimentset: str = attr.ib(kw_only=True, default=None)
    scoresets: List[str] = attr.ib(kw_only=True, default=None)

    def api_url() -> str:
        """

        Returns
        -------

        """
        return 'experiments/'

    def api_id_field() -> str:
        """

        Returns
        -------

        """
        return 'urn'

    def deserialize(json_dict):
        """
        Takes a json dictionary and returns an instance of this class.
        """
        return Experiment(**json_dict)


@attr.s
class NewExperiment(NewDataset):
    experimentset: str = attr.ib(kw_only=True, default=None)

    def api_url() -> str:
        """

        Returns
        -------

        """
        return 'experiments'

    def post_payload(self):
        """
        Use this to POST an instance of this class.
        """
        json_dict, files = prepare_for_encoding(
            attr.asdict(
                self,
                filter=attrs_filter,
                retain_collection_types=True,
                value_serializer=attrs_serializer
            )
        )
        return json_dict, files
