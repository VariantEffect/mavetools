import attr


@attr.s
class Licence:
    """
    Class contains attributes of Licence
    """

    long_name: str = attr.ib(kw_only=True, default=None)
    short_name: str = attr.ib(kw_only=True)
    link: str = attr.ib(kw_only=True, default=None)
    version: str = attr.ib(kw_only=True, default=None)

    def valid_short_names():
        """
        Returns list of valid Licence short names
        """
        return ["CC0", "CC BY-NC-SA 4.0", "CC BY 4.0"]
