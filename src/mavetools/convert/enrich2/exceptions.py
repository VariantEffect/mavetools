class HGVSValidationError(Exception):
    """
    Throw exception when a variant defines a SNP with a reference base
    that does not match that in its wild-type sequence/codon list.
    """

    pass


class HGVSMatchError(Exception):
    """
    Throw exception when a variant could not be pattern matched.
    """

    pass


class SequenceIndexError(Exception):
    """
    Throw exception when a variant defines a SNP with a position that goes
    out of bounds relative to its wild-type sequence.
    """

    pass


class InvalidWildTypeSequence(Exception):
    """
    Thrown if a wild-type sequence contains unknown characters.
    """


class SequenceFrameError(Exception):
    """
    Thrown if a sequence is not a multiple of three.
    """

    pass


class InvalidVariantType(Exception):
    """
    Throw exception when a specific type of event is expected (sub, del, etc)
    but not found.
    """
