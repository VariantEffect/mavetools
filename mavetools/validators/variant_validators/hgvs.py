from functools import partial
from typing import Optional, Union
import re

from mavehgvs import Variant, MaveHgvsParseError
from mavetools.validators.exceptions import ValidationError

# from core.utilities import is_null
# Used in CSV formatting
NA_value = "NA"
null_values_re = re.compile(
    r"\s+|none|nan|na|undefined|n/a|null|nil|{}".format(NA_value), flags=re.IGNORECASE
)


def is_null(value):
    """Returns True if a stripped/lowercase value in in `nan_col_values`."""
    value = str(value).strip().lower()
    return null_values_re.fullmatch(value) or not value


from mavetools.validators.for_variant_validators.constants import (
    hgvs_nt_column,
    hgvs_splice_column,
    hgvs_pro_column,
)

# from core.utilities import is_null
def is_null(value):
    """Returns True if a stripped/lowercase value in in `nan_col_values`."""
    value = str(value).strip().lower()
    return null_values_re.fullmatch(value) or not value


def validate_hgvs_string(
    value: Union[str, bytes],
    column: Optional[str] = None,
    splice_present: bool = False,
    targetseq: Optional[str] = None,
    relaxed_ordering: bool = False,
) -> Optional[str]:
    if is_null(value):
        return None

    if hasattr(value, "decode"):
        value = value.decode()
    if not isinstance(value, str):
        raise ValidationError(
            "Variant HGVS values input must be strings. "
            "'{}' has the type '{}'.".format(value, type(value).__name__)
        )

    if value.lower() == "_sy":
        raise ValidationError(
            "_sy is no longer supported and should be replaced by p.(=)"
        )
    elif value.lower() == "_wt":
        raise ValidationError(
            "_wt is no longer supported and should be replaced by (cgnp).="
        )

    try:
        variant = Variant(
            s=value, targetseq=targetseq, relaxed_ordering=relaxed_ordering
        )
    except MaveHgvsParseError as error:
        raise ValidationError(f"{value}: {str(error)}")

    prefix = variant.prefix.lower()
    if column in ("nt", hgvs_nt_column):
        if splice_present:
            if prefix not in "g":
                raise ValidationError(
                    f"'{value}' is not a genomic variant (prefix 'g.'). "
                    f"Nucleotide variants must be genomic if transcript "
                    f"variants are also defined."
                )
        else:
            if prefix not in "cn":
                raise ValidationError(
                    f"'{value}' is not a transcript variant. The accepted "
                    f"transcript variant prefixes are 'c.', 'n.'."
                )
    elif column in ("splice", hgvs_splice_column):
        if prefix not in "cn":
            raise ValidationError(
                f"'{value}' is not a transcript variant. The accepted "
                f"transcript variant prefixes are 'c.', 'n.'."
            )
    elif column in ("p", hgvs_pro_column):
        if prefix not in "p":
            raise ValidationError(
                f"'{value}' is not a protein variant. The accepted "
                f"protein variant prefix is 'p.'."
            )
    else:
        raise ValueError("Unknown column '{}'. Expected nt, splice or p".format(column))

    return str(variant)


validate_nt_variant = partial(validate_hgvs_string, **{"column": "nt"})
validate_splice_variant = partial(validate_hgvs_string, **{"column": "splice"})
validate_pro_variant = partial(validate_hgvs_string, **{"column": "p"})
