import attr
from typing import Any, Dict, List, Optional

from .external_identifier import ExternalIdentifier


@attr.s
class Urn:
    """
    Instantiates the Urn object and declares urn attribute
    """

    urn: str = attr.ib(kw_only=True)


@attr.s
class Keyword:
    """
    Instantiates the Keyword object and declares text attribute
    """

    text: str = attr.ib(kw_only=True)


@attr.s
class TimeStamped:
    """
    Instantiates TimeStamped object and declares attributes
    """

    creation_date: str = attr.ib(kw_only=True)
    modification_date: str = attr.ib(kw_only=True)


@attr.s
class Dataset(TimeStamped, Urn):
    """
    Instantiates Dataset object and declares attributes
    """

    publish_date: str = attr.ib(kw_only=True)
    created_by: str = attr.ib(kw_only=True)
    modified_by: str = attr.ib(kw_only=True)
    # optional attributes
    approved: Optional[str] = attr.ib(kw_only=True, default=None)
    private: Optional[bool] = attr.ib(kw_only=True, default=None)
    last_child_value: Optional[Any] = attr.ib(kw_only=True, default=None)
    # optional attribute
    extra_metadata: Optional[Dict[str, str]] = attr.ib(kw_only=True, default=None)
    abstract_text: str = attr.ib(kw_only=True)
    method_text: str = attr.ib(kw_only=True)
    short_description: str = attr.ib(kw_only=True)
    title: str = attr.ib(kw_only=True)
    keywords: List[Keyword] = attr.ib(kw_only=True)
    sra_ids: Optional[List[ExternalIdentifier]] = attr.ib(kw_only=True, default=None)
    doi_ids: Optional[List[ExternalIdentifier]] = attr.ib(kw_only=True, default=None)
    pubmed_ids: Optional[List[ExternalIdentifier]] = attr.ib(kw_only=True, default=None)
    contributors: List[str] = attr.ib(kw_only=True)

    def deserialize():
        # why pass?
        pass


@attr.s
class NewDataset:
    """
    Instantiates NewDataset object and declares object attributes
    """

    title: str = attr.ib(kw_only=True)
    short_description: str = attr.ib(kw_only=True)
    # optional attributes
    abstract_text: Optional[str] = attr.ib(kw_only=True, default=None)
    method_text: Optional[str] = attr.ib(kw_only=True, default=None)
    keywords: Optional[List[Keyword]] = attr.ib(kw_only=True, default=None)
    # TODO: change this once you know what this is supposed to be
    doi_ids: Optional[List[str]] = attr.ib(kw_only=True, default=None)
    sra_ids: Optional[List[str]] = attr.ib(kw_only=True, default=None)
    pubmed_ids: Optional[List[str]] = attr.ib(kw_only=True, default=None)
