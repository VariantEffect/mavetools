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
class SequenceIdentifier():
    """
    Instantiates SequenceIdentifier and declares attributes
    """

    identifier: ExternalIdentifier = attr.ib(kw_only=True)
    offset: int = attr.ib(kw_only=True)


@attr.s
class Taxonomy:
    """
    Instantiates Taxonomy and declares attributes
    """

    taxId: str = attr.ib(kw_only=True)
    organismName: str = attr.ib(kw_only=True, default=None)
    commonName: str = attr.ib(kw_only=True)
    rank: str = attr.ib(kw_only=True)
    hasDescribedSpeciesName: bool = attr.ib(kw_only=True)
    articleReference: str = attr.ib(kw_only=True)
    genomeId: str = attr.ib(kw_only=True)
    id: int = attr.ib(kw_only=True)
    url: str = attr.ib(kw_only=True)

@attr.s
class TargetSequence:
    """
    Instantiates targetSequence object and declares attributes
    """

    sequenceType: str = attr.ib(kw_only=True)
    sequence: str = attr.ib(kw_only=True)
    label: str = attr.ib(kw_only=True)
    taxonomy: Taxonomy = attr.ib(kw_only=True)


@attr.s
class Target:
    """
    Instantiates target object and declares attributes
    """

    name: str = attr.ib(kw_only=True)
    category: str = attr.ib(kw_only=True)
    ExternalIdentifiers: List[SequenceIdentifier] = attr.ib(kw_only=True)
    reference_sequence: ReferenceSequence = attr.ib(kw_only=True)
    id: int = attr.ib(kw_only=True)
    targetSequence: TargetSequence = attr.ib(kw_only=True)
    targetAccession: str = attr.ib(kw_only=True)


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
