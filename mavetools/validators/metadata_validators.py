import re
import idutils

from mavetools.validators.exceptions import ValidationError

NA_value = "NA"
null_values_re = re.compile(
    r"\s+|none|nan|na|undefined|n/a|null|nil|{}".format(NA_value), flags=re.IGNORECASE
)


def is_null(value):
    """Returns True if a stripped/lowercase value in in `nan_col_values`."""
    value = str(value).strip().lower()
    return null_values_re.fullmatch(value) or not value


def validate_sra_identifier(identifier):
    from mavetools.validators.for_metadata_validators.metadata_models import (
        ExternalIdentifier,
    )

    if isinstance(identifier, ExternalIdentifier):
        identifier = identifier.identifier
    if not (
        idutils.is_sra(identifier)
        or idutils.is_bioproject(identifier)
        or idutils.is_geo(identifier)
        or idutils.is_arrayexpress_array(identifier)
        or idutils.is_arrayexpress_experiment(identifier)
    ):
        raise ValidationError(
            (
                "%(id)s is not a valid SRA, GEO, ArrayExpress or BioProject "
                "accession."
            ),
            params={"id": identifier},
        )


def validate_keyword(kw):
    from mavetools.validators.for_metadata_validators.metadata_models import Keyword

    if isinstance(kw, Keyword):
        kw = kw.text
    if is_null(kw) or not isinstance(kw, str):
        raise ValidationError(
            "%(kw)s is not a valid keyword. Keywords must be valid strings.",
            params={"kw": kw},
        )


def validate_pubmed_identifier(identifier):
    from mavetools.validators.for_metadata_validators.metadata_models import (
        ExternalIdentifier,
    )

    if isinstance(identifier, ExternalIdentifier):
        identifier = identifier.identifier
    if not idutils.is_pmid(identifier):
        raise ValidationError(
            "%(id)s is not a valid PubMed identifier.", params={"id": identifier}
        )


def validate_doi_identifier(identifier):
    from mavetools.validators.for_metadata_validators.metadata_models import (
        ExternalIdentifier,
    )

    print(type(identifier))
    if isinstance(identifier, ExternalIdentifier):
        identifier = identifier.identifier
    if not idutils.is_doi(identifier):
        raise ValidationError("%(id)s is not a valid DOI.", params={"id": identifier})


def validate_ensembl_identifier(identifier):
    from mavetools.validators.for_metadata_validators.metadata_models import (
        ExternalIdentifier,
    )

    if isinstance(identifier, ExternalIdentifier):
        identifier = identifier.identifier
    if not idutils.is_ensembl(identifier):
        raise ValidationError(
            "%(id)s is not a valid Ensembl accession.", params={"id": identifier}
        )


def validate_uniprot_identifier(identifier):
    from mavetools.validators.for_metadata_validators.metadata_models import (
        ExternalIdentifier,
    )

    if isinstance(identifier, ExternalIdentifier):
        identifier = identifier.identifier
    if not idutils.is_uniprot(identifier):
        raise ValidationError(
            "%(id)s is not a valid UniProt accession.", params={"id": identifier}
        )


def validate_refseq_identifier(identifier):
    from mavetools.validators.for_metadata_validators.metadata_models import (
        ExternalIdentifier,
    )

    if isinstance(identifier, ExternalIdentifier):
        identifier = identifier.identifier
    if not idutils.is_refseq(identifier):
        raise ValidationError(
            "%(id)s is not a valid RefSeq accession.", params={"id": identifier}
        )


def validate_genome_identifier(identifier):
    from mavetools.validators.for_metadata_validators.metadata_models import (
        ExternalIdentifier,
    )

    if isinstance(identifier, ExternalIdentifier):
        identifier = identifier.identifier
    if not idutils.is_genome(identifier):
        raise ValidationError(
            "%(id)s is not a valid GenBank or RefSeq genome assembly.",
            params={"id": identifier},
        )


def validate_keyword_list(values):
    for value in values:
        if not is_null(value):
            validate_keyword(value)


def validate_pubmed_list(values):
    for value in values:
        if not is_null(value):
            validate_pubmed_identifier(value)


def validate_sra_list(values):
    for value in values:
        if not is_null(value):
            validate_sra_identifier(value)


def validate_doi_list(values):
    for value in values:
        if not is_null(value):
            validate_doi_identifier(value)


def validate_ensembl_list(values):
    for value in values:
        if not is_null(value):
            validate_ensembl_identifier(value)


def validate_refseq_list(values):
    for value in values:
        if not is_null(value):
            validate_refseq_identifier(value)


def validate_uniprot_list(values):
    for value in values:
        if not is_null(value):
            validate_uniprot_identifier(value)
