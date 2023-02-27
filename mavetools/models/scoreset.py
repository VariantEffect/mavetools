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
    """
    This class instantiates the ScoreSet object and declares the fields of that object.
    It inherits the attributes of APIObject and Dataset
    """

    experiment: str = attr.ib(kw_only=True)
    licenceId: Licence = attr.ib(kw_only=True)
    targetGene: Target = attr.ib(kw_only=True)
    # optional attributes
    datasetColumns: Optional[Any] = attr.ib(kw_only=True, default=None)
    replacesId: Optional[str] = attr.ib(kw_only=True, default=None)
    #score_columns: Optional[List[str] = attr.ib(kw_only=True)#TODO see how these should be structured within DatasetColumns
    #count_columns: Optional[List[str]] = attr.ib(kw_only=True)#
    # optional attributes
    previous_version: Optional[str] = attr.ib(kw_only=True, default=None)
    next_version: Optional[str] = attr.ib(kw_only=True, default=None)
    #current_version: Optional[str] = attr.ib(kw_only=True)
    numVariants: int = attr.ib(kw_only=True)
    dataUsagePolicy: str = attr.ib(kw_only=True)
    is_meta_analysis: Optional[bool] = attr.ib(kw_only=True, default=None)

    def api_url() -> str:
        """
        Returns API endpoint
        """
        return "scoresets/"

    def api_id_field() -> str:
        """
        Returns API ID field as string urn
        """
        return "urn"

    def deserialize(json_dict):
        """
        Takes a json dictionary and returns an instance of this class.
        """
        return ScoreSet(**json_dict)


@attr.s
class NewScoreSet(NewDataset):
    """
    This Class instantiates NewScoreSet and declares the variables present in NewScoreSet object.
    Attributes are set before posting a model instance
    """

    experiment: str = attr.ib(kw_only=True)

    # the following fields are optional
    meta_analysis_for: Optional[str] = attr.ib(kw_only=True, default=None)
    replacesId: Optional[str] = attr.ib(kw_only=True, default=None)
    licenceId: Optional[str] = attr.ib(kw_only=True, default=None)
    dataUsagePolicy: Optional[str] = attr.ib(kw_only=True, default=None)

    # These can be strings or open filepaths
    score_data: Union[str, BinaryIO] = attr.ib(kw_only=True)
    count_data: Union[str, BinaryIO] = attr.ib(kw_only=True)
    meta_data: Union[str, BinaryIO] = attr.ib(kw_only=True)


@attr.s
class NewScoreSetRequest(APIObject):
    """
    This class instantiates NewScoreSetRequest and sets the fields for NewScoreSetRequest.
    Attributes are set before posting a model instance
    """

    scoreset: NewScoreSet = attr.ib(kw_only=True)
    targetGene: NewTarget = attr.ib(kw_only=True)

    # These only need the ExternalIdentifier identifier field
    uniprot: SequenceOffset = attr.ib(kw_only=True, default=None)
    ensembl: SequenceOffset = attr.ib(kw_only=True, default=None)
    refseq: SequenceOffset = attr.ib(kw_only=True, default=None)

    reference_maps: List[ReferenceMap] = attr.ib(kw_only=True)

    def api_url() -> str:
        """
        Returns API endpoint
        """
        return "scoresets"

    def post_payload(self):
        """
        Use this to POST an instance of this class.
        data is converted to appropriate type and returned
        """
        json_dict, files = prepare_for_encoding(
            attr.asdict(
                self,
                filter=attrs_filter,
                retain_collection_types=True,
                value_serializer=attrs_serializer,
            )
        )
        return json_dict, files
