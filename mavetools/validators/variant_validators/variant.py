from typing import Dict

from mavetools.validators.for_variant_validators.constants import (
    variant_score_data,
    variant_count_data,
    required_score_column,
)
from mavetools.validators.exceptions import ValidationError


def validate_columns_match(variant, scoreset) -> None:
    """
    Validate that a child matches parents defined columns to keep
    data in sync.
    """
    try:
        if variant.score_columns != scoreset.score_columns:
            raise ValidationError(
                f"Variant defines score columns '{variant.score_columns}' "
                f"but parent defines columns '{scoreset.score_columns}. "
            )
        if variant.count_columns != scoreset.count_columns:
            raise ValidationError(
                f"Variant defines count columns '{variant.count_columns}' "
                f"but parent defines columns '{scoreset.count_columns}. "
            )
    except KeyError as error:
        raise ValidationError(f"Missing key {str(error)}")


def validate_variant_json(data: Dict[str, Dict]) -> None:
    """
    Checks a given dictionary to ensure that it is suitable to be used
    as the `data` attribute in a :class:`Variant` instance.

    Parameters
    ----------
    data : dict
        Dictionary of keys mapping to a list.
    """
    expected_keys = [variant_score_data, variant_count_data]
    for key in expected_keys:
        if key not in data.keys():
            raise ValidationError(f"Missing the required key {key}")

    if required_score_column not in data[variant_score_data]:
        raise ValidationError(
            f"Missing required column '{required_score_column}' in variant's score data."
        )

    extras = [k for k in data.keys() if k not in set(expected_keys)]
    if len(extras) > 0:
        extras = [k for k in data.keys() if k not in expected_keys]
        raise ValidationError("Encountered unexpected keys {extras}")

    # Check the correct data types are given.
    for key in expected_keys:
        if not isinstance(data[key], dict):
            type_ = type(data[key]).__name__
            raise ValidationError(f"Value for '{key}' must be a dict not {type_}.")
