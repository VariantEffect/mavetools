import attr
import json
from typing import List, Optional

from .external_identifier import ExternalIdentifier


@attr.s
class ReferenceSequence:
    """
    Instantiates ReferenceSequence object and declares attributes
    """

    sequence: str = attr.ib(kw_only=True)


@attr.s
class SequenceOffset(ExternalIdentifier):
    """
    Instantiates SequenceOffset and declares attributes
    """

    offset: int = attr.ib(kw_only=True)


@attr.s
class ReferenceGenome:
    """
    Instantiates SequenceOffset and declares attributes
    """

    short_name: str = attr.ib(kw_only=True)
    organism_name: str = attr.ib(kw_only=True, default=None)
    assembly_identifier: ExternalIdentifier = attr.ib(kw_only=True, default=None)


@attr.s
class ReferenceMap:
    """
    Instantiates ReferenceMap object and declares genome attribute
    """

    genome: ReferenceGenome = attr.ib(kw_only=True)


@attr.s
class Target:
    """
    Instantiates target object and declares attributes
    """

    name: str = attr.ib(kw_only=True)
    reference_sequence: ReferenceSequence = attr.ib(kw_only=True)
    uniprot: SequenceOffset = attr.ib(kw_only=True)
    ensembl: SequenceOffset = attr.ib(kw_only=True)
    refseq: SequenceOffset = attr.ib(kw_only=True)

    reference_maps: List[ReferenceMap] = attr.ib(kw_only=True)
    scoreset: str = attr.ib(kw_only=True)
    type: str = attr.ib(kw_only=True)


@attr.s
class NewTarget:
    """
    Instantiates NewTarget object and declares attributes and functions
    """

    # existing_target: str = attr.ib(kw_only=True, default=None)
    name: str = attr.ib(kw_only=True)
    reference_sequence: Optional[ReferenceSequence] = attr.ib(
        kw_only=True, default=None
    )
    type: str = attr.ib(kw_only=True)
    sequence_type: str = attr.ib(kw_only=True)
    sequence_text: str = attr.ib(kw_only=True, default=None)

    # Needs to be a file path
    fasta_file: str = attr.ib(kw_only=True, default=None)

    def valid_types():
        """
        Returns list of valid target types
        """
        return ["Protein coding", "Regulatory", "Other noncoding"]

    def valid_sequence_types():
        """
        Returns list of valid sequence types
        """
        return ["Infer", "DNA", "Protein"]
