from mavetools.validators import dataset_validators
from mavetools.validators import urn_validators


def validate_all(countfile=None, scorefile=None, scorejson=None,
                 scoreset_urn=None, experiment_urn=None, experimentset_urn=None):
    """
    By calling other helper functions, this function runs all of the validation code
    """
    validate_dataset(countfile, scorefile, scorejson)
    validate_urn(scoreset_urn, experiment_urn, experimentset_urn)


def validate_dataset(countfile=None, scorefile=None, scorejson=None):
    """
    This function calls all of the validation functions within
    mavetools/mavetools/validators/dataset_validation.py

    Parameters
    ----------
    countfile
    scorefile
    scorejson

    Returns
    -------


    """

    # how to incorporate word limit validator?

    if scorefile is not None:
        # open scorefile
        file_scorefile = open(scorefile, "r")

        # this one returns header
        #scoreheader = dataset_validators.read_header_from_io(file=scorefile)

        # if the header was returned, do these ones
        #dataset_validators.validate_has_hgvs_in_header(scorefile)
        #dataset_validators.validate_at_least_one_additional_column(header=scoreheader)
        #dataset_validators.validate_header_contains_no_null_columns(header=scoreheader)

        dataset_validators.validate_scoreset_score_data_input(file=file_scorefile)

        if scorejson is not None:
            # open scorejson
            file_scorejson = open(scorejson, "r")
            dataset_validators.validate_scoreset_json(dict_=file_scorejson)

    if countfile is not None:
        # open countfile
        file_countfile = open(countfile, "r")
        countheader = dataset_validators.read_header_from_io(file=file_countfile)

        # if the header was returned, do these ones
        #dataset_validators.validate_has_hgvs_in_header(header=countheader)
        #dataset_validators.validate_at_least_one_additional_column(header=countheader)
        #dataset_validators.validate_header_contains_no_null_columns(header=countheader)

        dataset_validators.validate_scoreset_count_data_input(file=file_countfile)

    #if scorefile is not None and countfile is not None:
        #dataset_validators.validate_datasets_define_same_variants(scores=file_scorefile, counts=file_countfile)


def validate_urn(scoreset_urn=None, experiment_urn=None, experimentset_urn=None):
    """

    Parameters
    ----------
    scoreset_urn
    experiment_urn
    experimentset_urn

    Returns
    -------

    """
    if scoreset_urn is not None:
        urn_validators.validate_mavedb_urn(scoreset_urn)
        #urn_validators.validate_mavedb_urn_variant(scoreset_urn)
        urn_validators.validate_mavedb_urn_scoreset(scoreset_urn)
    if experiment_urn is not None:
        urn_validators.validate_mavedb_urn(experiment_urn)
        #urn_validators.validate_mavedb_urn_variant(experiment_urn)
        urn_validators.validate_mavedb_urn_experiment(experiment_urn)
    if experimentset_urn is not None:
        urn_validators.validate_mavedb_urn(experimentset_urn)
        #urn_validators.validate_mavedb_urn_variant(experimentset_urn)
        urn_validators.validate_mavedb_urn_experimentset(experimentset_urn)