import datetime
from typing import Union, Optional

from django.contrib.auth import get_user_model
from django.db import transaction

from dataset import models
from variant.models import Variant
from urn.models import get_model_by_urn

from .models.experimentset import ExperimentSet
from .models.experiment import Experiment
from .models.scoreset import ScoreSet

User = get_user_model()


@transaction.atomic
def publish_dataset(
    dataset: Union[ExperimentSet, Experiment, ScoreSet], user: Optional[User] = None
) -> Union[ExperimentSet, Experiment, ScoreSet]:
    """
    Publishes a dataset by traversing the parent tree. Assigns a public
    urn of the format <urn:mavedb:X>, sets the private bit and associated
    publish metadata.

    Does nothing if the dataset already has a public urn or is not private.

    Parameters
    ----------
    dataset : ExperimentSet | Experiment | ScoreSet
        The dataset to publish.

    user : User
        The user requesting the publish.

    Raises
    ------
    TypeError : Not a dataset
    """
    if not isinstance(dataset, (ExperimentSet, Experiment, ScoreSet)):
        raise TypeError(
            "Expected a ExperimentSet, Experiment or ScoreSet instance. "
            "Found {}".format(dataset.__class__.__name__)
        )

    if not dataset.private or dataset.has_public_urn:
        return dataset

    scoreset = None
    experiment = None
    # Forces a full refresh on on the dataset including nested parents.
    dataset = get_model_by_urn(dataset.urn)

    if isinstance(dataset, models.scoreset.ScoreSet):
        experimentset = models.experimentset.assign_public_urn(
            dataset.experiment.experimentset
        )
        experiment = models.experiment.assign_public_urn(dataset.experiment)
        scoreset = models.scoreset.assign_public_urn(dataset)
        urns = Variant.bulk_create_urns(
            scoreset.children.count(), scoreset, reset_counter=True
        )
        for urn, child in zip(urns, scoreset.children.all()):
            child.urn = urn
            child.save()
    elif isinstance(dataset, models.experiment.Experiment):
        experimentset = models.experimentset.assign_public_urn(dataset.experimentset)
        experiment = models.experiment.assign_public_urn(dataset)
    elif isinstance(dataset, models.experimentset.ExperimentSet):
        experimentset = models.experimentset.assign_public_urn(dataset)
    else:
        raise TypeError(
            "Expected ExperimentSet, Experiment or ScoreSet. Found {}".format(
                dataset.__class__.__name__
            )
        )
    # Note: assigning a public urn to a child will alter `last_child_value`
    # of the parent. Before saving below, call `refresh_from_db` to update
    # any changes made when a child urn is assigned. Otherwise an outdated
    # version of the parent will over-write the most recent changes.
    if scoreset:
        scoreset.refresh_from_db()
        scoreset.publish_date = datetime.date.today()
        scoreset.private = False
        scoreset.set_modified_by(user, propagate=False)
        scoreset.save()

    if experiment:
        experiment.refresh_from_db()
        experiment.publish_date = datetime.date.today()
        experiment.private = False
        experiment.set_modified_by(user, propagate=False)
        experiment.save()

    if experimentset:
        experimentset.refresh_from_db()
        experimentset.publish_date = datetime.date.today()
        experimentset.private = False
        experimentset.set_modified_by(user, propagate=False)
        experimentset.save()

    dataset.refresh_from_db()
    return get_model_by_urn(dataset.urn)  # Full refresh on nested parents.
