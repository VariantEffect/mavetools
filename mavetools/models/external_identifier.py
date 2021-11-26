import attr
from typing import Optional

@attr.s
class ExternalIdentifier():
    identifier: str = attr.ib(kw_only=True)
    url: Optional[str] = attr.ib(kw_only=True, default=None)
    dbversion: Optional[str] = attr.ib(kw_only=True, default=None)
    dbname: Optional[str] = attr.ib(kw_only=True, default=None)
