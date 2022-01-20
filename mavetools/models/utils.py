import attr
import collections
import os
from typing import BinaryIO

from .licence import Licence


def attrs_filter(attr, value):
    """
    ???
    Parameters
    ----------
    attr
    value

    Returns
    -------

    """
    return value is not None


def attrs_serializer(inst, field, value):
    """
    ???
    Parameters
    ----------
    inst
    field
    value

    Returns
    -------

    """
    if isinstance(value, str):
        if os.path.isfile(value):
            ext = os.path.splitext(value)[1]
            return (f"{field.name}{ext}", open(value, "rb"), "application/octet-stream")
        return value
    if value is not None:
        return value


def prepare_for_encoding(nested_dict):
    """
    Prepares data for encoding by converting the data in the provided nested_dict into
    a json_dict and file_dict
    Parameters
    ----------
    nested_dict (dictionary): data to be converted

    Returns
    -------
    json_dict
    file_dict
    """
    json_dict = {}
    file_dict = {}
    for k, v in nested_dict.items():
        if isinstance(v, tuple):
            file_dict[k] = v
        elif isinstance(v, collections.MutableMapping):
            j, f = prepare_for_encoding(v)
            # Keep the original keys for nested dicts, but flatten the files dict
            json_dict[k] = j
            file_dict.update(f)
        else:
            json_dict[k] = v

    return json_dict, file_dict
