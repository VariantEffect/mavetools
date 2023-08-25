from typing import Mapping, Optional
from mavedb.view_models.score_set import ScoreSetCreate
from mavedb.view_models.experiment import ExperimentCreate
import humps


def infer_record_type(record: Mapping) -> Optional[str]:
    """
    Infer the type of the record.

    Parameters
    ----------
    record : Mapping
        The record to infer as a dictionary or other mapping type.

    Returns
    -------
    Optional[str]
        One of "experiment" or "score_set" if the dataset type can be inferred; else None.
    """
    # TODO: make this more specific
    if all(x in humps.decamelize(record) for x in ("title", "target_gene")):
        return "score_set"
    elif "title" in record.keys():
        return "experiment"
    else:
        return None


def validate_dataset_with_create_model(dataset: Mapping) -> None:
    """
    Validate a dataset using a MaveDB view model before uploading it to the API.

    Parameters
    ----------
    dataset : Mapping
        The dataset in a dictionary-like format.

    Returns
    -------
    None

    Raises
    ------
    ValueError
        If no suitable MaveDB view model was available.

    """
    record_type = infer_record_type(dataset)
    if record_type == "score_set":
        ScoreSetCreate(**dataset)
    elif record_type == "experiment":
        ExperimentCreate(**dataset)
    else:
        raise ValueError("could not find a model for validating the dataset")
