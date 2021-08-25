import attr
import json
import os
from typing import Any, BinaryIO, Dict, List, Optional, Union

from .base import APIObject
from .dataset import Dataset, NewDataset
from .licence import Licence
from .target import NewTarget, ReferenceMap, SequenceOffset, Target
from .utils import attrs_filter, attrs_serializer, prepare_for_encoding


@attr.s
class ScoreSet(APIObject, Dataset):
    experiment: str = attr.ib(kw_only=True)
    licence: Licence = attr.ib(kw_only=True)
    target: Target = attr.ib(kw_only=True)
    dataset_columns: Optional[Any] = attr.ib(kw_only=True, default=None)
    replaces: Optional[str] = attr.ib(kw_only=True, default=None)
    score_columns: List[str] = attr.ib(kw_only=True)
    count_columns: List[str] = attr.ib(kw_only=True)
    previous_version: Optional[str] = attr.ib(kw_only=True, default=None)
    next_version: Optional[str] = attr.ib(kw_only=True, default=None)
    current_version: str = attr.ib(kw_only=True)
    variant_count: int = attr.ib(kw_only=True)
    data_usage_policy: str = attr.ib(kw_only=True)
    is_meta_analysis: bool = attr.ib(kw_only=True)

    def api_url() -> str:
        """

        Returns
        -------
        API endpoint
        """
        return 'scoresets/'

    def api_id_field() -> str:
        return 'urn'

    def deserialize(json_dict):
        """
        Takes a json dictionary and returns an instance of this class.
        """
        return ScoreSet(**json_dict)


@attr.s
class NewScoreSet(NewDataset):
    experiment: str = attr.ib(kw_only=True)
    meta_analysis_for: Optional[str] = attr.ib(kw_only=True, default=None)
    replaces: Optional[str] = attr.ib(kw_only=True, default=None)
    licence: Optional[str] = attr.ib(kw_only=True, default=None)
    data_usage_policy: Optional[str] = attr.ib(kw_only=True, default=None)

    # These can be strings or open filepaths
    score_data: Union[str, BinaryIO] = attr.ib(kw_only=True)
    count_data: Union[str, BinaryIO] = attr.ib(kw_only=True)
    meta_data: Union[str, BinaryIO] = attr.ib(kw_only=True)


@attr.s
class NewScoreSetRequest(APIObject):
    scoreset: NewScoreSet = attr.ib(kw_only=True)
    target: NewTarget = attr.ib(kw_only=True)

    # These only need the ExternalIdentifier identifier field
    uniprot: SequenceOffset = attr.ib(kw_only=True, default=None)
    ensembl: SequenceOffset = attr.ib(kw_only=True, default=None)
    refseq: SequenceOffset = attr.ib(kw_only=True, default=None)

    reference_maps: List[ReferenceMap] = attr.ib(kw_only=True)

    def api_url() -> str:
        return 'scoresets'

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
