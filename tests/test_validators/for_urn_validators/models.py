import random
import string

from mavetools.validators.urn_validators import (
    MAVEDB_EXPERIMENTSET_URN_DIGITS,
    MAVEDB_URN_MAX_LENGTH,
    MAVEDB_URN_NAMESPACE,
    MAVEDB_TMP_URN_DIGITS,
)


RANDOM_CHARS = string.ascii_lowercase + string.ascii_uppercase + string.digits


def generate_tmp_urn():
    return "tmp:{}".format(
        "".join([random.choice(RANDOM_CHARS) for _ in range(MAVEDB_TMP_URN_DIGITS)])
    )
