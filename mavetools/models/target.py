import attr
import json
from typing import List, Optional

from .external_identifier import ExternalIdentifier


@attr.s
class ReferenceSequence():
    sequence: str = attr.ib(kw_only=True)


@attr.s
class SequenceOffset(ExternalIdentifier):
    offset: int = attr.ib(kw_only=True)


@attr.s
class ReferenceGenome():
    short_name: str = attr.ib(kw_only=True)
    organism_name: str = attr.ib(kw_only=True, default=None)
    assembly_identifier: ExternalIdentifier = attr.ib(kw_only=True, default=None)


@attr.s
class ReferenceMap():
    genome: ReferenceGenome = attr.ib(kw_only=True)


@attr.s
class Target():
    name: str = attr.ib(kw_only=True)
    reference_sequence: ReferenceSequence = attr.ib(kw_only=True)
    uniprot: SequenceOffset = attr.ib(kw_only=True)
    ensembl: SequenceOffset = attr.ib(kw_only=True)
    refseq: SequenceOffset = attr.ib(kw_only=True)

    reference_maps: List[ReferenceMap] = attr.ib(kw_only=True)
    scoreset: str = attr.ib(kw_only=True)
    type: str = attr.ib(kw_only=True)


@attr.s
class NewTarget():
    # existing_target: str = attr.ib(kw_only=True, default=None)
    name: str = attr.ib(kw_only=True)
    reference_sequence: Optional[ReferenceSequence] = attr.ib(kw_only=True, default=None)
    type: str = attr.ib(kw_only=True)
    sequence_type: str = attr.ib(kw_only=True)
    sequence_text: str = attr.ib(kw_only=True, default=None)

    # Needs to be a file path
    fasta_file: str = attr.ib(kw_only=True, default=None)

    def valid_types():
        return ['Protein coding', 'Regulatory', 'Other noncoding']

    def valid_sequence_types():
        return ['Infer', 'DNA', 'Protein']
