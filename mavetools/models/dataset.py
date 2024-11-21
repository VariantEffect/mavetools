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

    creationDate: str = attr.ib(kw_only=True)
    modificationDate: str = attr.ib(kw_only=True)


@attr.s
class Dataset(TimeStamped, Urn):
    """
    Instantiates Dataset object and declares attributes
    """
    # record keeping attributes
    publishedDate: str = attr.ib(kw_only=True)
    createdBy: str = attr.ib(kw_only=True)
    modifiedBy: str = attr.ib(kw_only=True)

    # required attributes
    shortDescription: str = attr.ib(kw_only=True)
    title: str = attr.ib(kw_only=True)


    # optional attribute
    contributors: Optional[List[str]] = attr.ib(kw_only=True, default=None)
    abstractText: str = attr.ib(kw_only=True)
    methodText: str = attr.ib(kw_only=True)
    keywords: List[Keyword] = attr.ib(kw_only=True)
    secondaryPublicationIdentifiers: Optional[List[ExternalIdentifier]] = attr.ib(kw_only=True, default=None)
    doiIdentifiers: Optional[List[ExternalIdentifier]] = attr.ib(kw_only=True, default=None)
    primaryPublicationIdentifiers: Optional[List[ExternalIdentifier]] = attr.ib(kw_only=True, default=None)
    extraMetadata: Optional[Dict[str, str]] = attr.ib(kw_only=True, default=None)
    private: Optional[bool] = attr.ib(kw_only=True, default=None)
    processingState: Optional[str] = attr.ib(kw_only=True, default=None)
    processingErrors: Optional[str] = attr.ib(kw_only=True, default=None)

    def deserialize():
        raise NotImplementedError()
