import attr
from typing import Optional


@attr.s
class ExternalIdentifier:
    """
    Instantiates ExternalIdentifier object and declares object attributes
    """

    identifier: str = attr.ib(kw_only=True)
    # optional attributes
    url: Optional[str] = attr.ib(kw_only=True, default=None)
    dbversion: Optional[str] = attr.ib(kw_only=True, default=None)
    dbname: Optional[str] = attr.ib(kw_only=True, default=None)
