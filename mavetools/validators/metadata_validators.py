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
    if not (
        idutils.is_sra(identifier)
        or idutils.is_bioproject(identifier)
        or idutils.is_geo(identifier)
        or idutils.is_arrayexpress_array(identifier)
        or idutils.is_arrayexpress_experiment(identifier)
    ):
        raise ValidationError(
            f"'{identifier} is not a valid SRA, GEO, ArrayExpress or BioProject "
            "accession."
        )


def validate_keyword(kw):
    if is_null(kw) or not isinstance(kw, str):
        raise ValidationError(
            f"'{kw}' not a valid keyword. Keywords must be valid strings."
        )


def validate_pubmed_identifier(identifier):
    if not idutils.is_pmid(identifier):
        raise ValidationError(f"'{identifier} is not a valid PubMed identifier.")


def validate_doi_identifier(identifier):
    if not idutils.is_doi(identifier):
        raise ValidationError(f"'{identifier}' is not a valid DOI.")


def validate_ensembl_identifier(identifier):
    if not idutils.is_ensembl(identifier):
        raise ValidationError(f"'{identifier}' is not a valid Ensembl accession.")


def validate_uniprot_identifier(identifier):
    if not idutils.is_uniprot(identifier):
        raise ValidationError(f"'{identifier}' is not a valid UniProt accession.")


def validate_refseq_identifier(identifier):
    if not idutils.is_refseq(identifier):
        raise ValidationError(f"'{identifier}' is not a valid RefSeq accession.")


def validate_genome_identifier(identifier):
    if not idutils.is_genome(identifier):
        raise ValidationError(
            f"'{identifier}' is not a valid GenBank or RefSeq genome assembly."
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
