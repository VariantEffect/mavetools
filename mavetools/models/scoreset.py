import attr
import json
import os
from typing import Any, BinaryIO, Dict, List, Optional, Union
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from .base import APIObject
from .dataset import Dataset, NewDataset
from .licence import Licence
from .target import NewTarget, ReferenceMap, SequenceOffset, Target
from .utils import attrs_filter, attrs_serializer, prepare_for_encoding, get_variant_type, median, score_scale_function
from fqfa.util.translate import translate_dna

from mavetools.models.sequence_retrieval import getUniprotSequence

@attr.s
class ScoreSet(APIObject, Dataset):
    """
    This class instantiates the ScoreSet object and declares the fields of that object.
    It inherits the attributes of APIObject and Dataset
    """

    experiment: str = attr.ib(kw_only=True)
    licence: Licence = attr.ib(kw_only=True)
    target: Target = attr.ib(kw_only=True)
    # optional attributes
    dataset_columns: Optional[Any] = attr.ib(kw_only=True, default=None)
    replaces: Optional[str] = attr.ib(kw_only=True, default=None)
    score_columns: List[str] = attr.ib(kw_only=True)
    count_columns: List[str] = attr.ib(kw_only=True)
    # optional attributes
    previous_version: Optional[str] = attr.ib(kw_only=True, default=None)
    next_version: Optional[str] = attr.ib(kw_only=True, default=None)
    current_version: str = attr.ib(kw_only=True)
    variant_count: int = attr.ib(kw_only=True)
    data_usage_policy: str = attr.ib(kw_only=True)
    is_meta_analysis: bool = attr.ib(kw_only=True)


    def api_url() -> str:
        """
        Returns API endpoint
        """
        return "scoresets/"

    def api_id_field() -> str:
        """
        Returns API ID field as string urn
        """
        return "urn"

    def deserialize(json_dict):
        """
        Takes a json dictionary and returns an instance of this class.
        """
        return ScoreSet(**json_dict)

    def get_protein_sequence(self):
        raw_seq = self.target['reference_sequence']['sequence']
        if self.target['reference_sequence']['sequence_type'] == 'protein':
            return raw_seq
        else:
            return translate_dna(raw_seq)[0]

    def get_full_sequence_info(self):
        if self.target['uniprot'] != None:
            full_seq = getUniprotSequence(self.target['uniprot']['identifier'])
            offset = self.target['uniprot']['offset']
        elif self.target['refseq'] != None:
            full_seq = get_refseq_sequences(self.target['refseq']['identifier'])[self.target['refseq']['identifier']]
            offset = self.target['refseq']['offset']
        elif self.target['ensembl'] != None:
            print(f'No Uniprot, no RefSeq, but Ensembl ID: {self.urn}')
        else:
            full_seq = self.get_protein_sequence()
            offset = 0
        return full_seq, offset

    def get_score_table_positions(self):
        for pos, score_table_header in enumerate(self.score_columns):
            if score_table_header == 'hgvs_pro':
                hgvs_pro_pos = pos + 1
            if score_table_header == 'hgvs_nt':
                hgvs_nt_pos = pos + 1
            if score_table_header == 'score':
                score_pos = pos + 1
        return hgvs_pro_pos, hgvs_nt_pos, score_pos

class ScoreSetData:

    """
    This class represents the score tables for a corresponding scoresheet.
    Can also represent score tables for a corresponding experiment,
    in that case the scoretables from all scoresheets belonging to the experiment
    got aggregated (see aggregate_scoresetdata function of MLdataset class in ml_tools.py).
    """

    def __init__(self, urn, csv_formatted_data, hgvs_pro_pos = None, score_pos = None, hgvs_nt_pos = None):

        """
        Initialization routine that builds a class instance from the content of a MaveDB scoreset table.

        Parameters
        ----------

        urn
            MaveDB urn identifier

        csv_formatted_data
            Content of a MaveDB scoreset table.

        hgvs_pro_pos
            Column number of the column that contains the hgvs_pro identifiers.

        score_pos
            Column number of the column that contains the experimental effect values.

        hgvs_pro_pos
            Column number of the column that contains the hgvs_nt identifiers.

        Returns
        -------

        class_instance
            Instance of ScoreSetData class.
        """

        #Initialize several class attributes
        self.urn = urn
        self.score_dict = {}
        self.nt_score_dict = {}
        self.sav_scores = {}
        self.synonymous_scores = {}
        self.mean_score = None
        self.min_score = None
        self.max_score = None
        self.mean_sav = None
        self.min_sav = None
        self.max_sav = None
        self.cardinality = 0
        self.sav_cardinality = 0
        self.sorted_sav_values = None

        #Parse the table

        lines = csv_formatted_data.split('\n')
        for line in lines:
            if line == '':
                continue
            if line[0] == '#':
                continue
            words = line.split(',')

            #In case their is a header row, parse the hgvs_pro, hgvs_nt, and effect column numbers
            if words[0][:4] != 'urn:':
                for pos, word in enumerate(words):
                    if word == 'hgvs_pro':
                        hgvs_pro_pos = pos
                    if word == 'score':
                        score_pos = pos
                    if word == 'hgvs_nt':
                        hgvs_nt_pos = pos
            else:
                if words[score_pos] == 'NA':
                    continue
                try:
                    score = float(words[score_pos])
                except:
                    #Debugging help, can be removed, when better sanity checks are implemented
                    print(score_pos)
                    print(words[score_pos])
                    print(lines[10:])
                    sys.exit()

                #Retrieve the column values and write them to their corresponding class attribute datastructures
                hgvs_pro = words[hgvs_pro_pos]
                hgvs_nt = words[hgvs_nt_pos]

                if hgvs_pro != 'NA':
                    self.score_dict[hgvs_pro] = score
                elif hgvs_nt != 'NA':
                    self.nt_score_dict[hgvs_nt] = score
                else:
                    continue

        #Get some basic information of the data

        for hgvs_pro in self.score_dict:
            score = self.score_dict[hgvs_pro]
            self.cardinality += 1

            #Dependant on the variant type, we need to do different things later
            #For now, store separately synonymous variants and single amino acid variations
            #There is perhaps the need to implement similar datastructures for indels and non-coding variations.
            variant_type = get_variant_type(hgvs_pro)

            if self.min_score is None:
                self.min_score = score
            elif score < self.min_score:
                self.min_score = score

            if self.max_score is None:
                self.max_score = score
            elif score > self.max_score:
                self.max_score = score

            if variant_type == 'synonymous':
                self.synonymous_scores[hgvs_pro] = score
            elif variant_type == 'sav':
                if self.min_sav is None:
                    self.min_sav = score
                elif score < self.min_sav:
                    self.min_sav = score

                if self.max_sav is None:
                    self.max_sav = score
                elif score > self.max_sav:
                    self.max_sav = score

                self.sav_scores[hgvs_pro] = score
                self.sav_cardinality += 1

        if self.cardinality > 0:
            self.mean_score = sum(self.score_dict.values()) / self.cardinality

        if self.sav_cardinality > 0:
            self.mean_sav = sum(self.sav_scores.values()) / self.sav_cardinality


    def get_sorted_sav_values(self):

        """
        Sorts the SAV effect values and stores them, so if called again, the sorting don't has to be recalculated again.

        Returns
        -------

        self.sorted_sav_values
            Pointer to the sorted sav value datastructure.
        """

        if self.sorted_sav_values is not None:
            return self.sorted_sav_values
        else:
            self.sorted_sav_values = sorted(self.sav_scores.values())
            return self.sorted_sav_values

    def sav_binning(self):

        """
        Creates a 10-bin histogram of all sav scores.

        Returns
        -------

        bins
            The calculated histogram
        """

        n_of_bins = 10

        sorted_scores = self.get_sorted_sav_values()
        bin_size = (self.max_sav - self.min_sav)/n_of_bins

        bins = []
        for i in range(n_of_bins):
            bins.append([])
        current_bin_number = 0
        for score in sorted_scores:
            if score > (self.min_sav + bin_size*(current_bin_number+1)):
                current_bin_number += 1
            try:
                bins[current_bin_number].append(score)
            except:
                bins[current_bin_number-1].append(score)

        return bins


    def scale_sav_data(self, verbosity = 0):

        """
        Scales all sav scores. Scaling is based on a typical neutral effect score and a typical strong effect value.
        For more details of the scaling process, look into the score_scale_function in utils.py.


        Parameters
        ---------

        verbosity
            Verbosity level of the function

        Returns
        -------

        None
            if the typical neutral effect and the strong effect value are equal, the scaling would be non-sensical
        """

        median_synonymous = median(self.synonymous_scores.values())
        sorted_scores = self.get_sorted_sav_values()
        perc_size = max([1,round(self.sav_cardinality/100)])

        self.scaled_sav_scores = {}

        if median_synonymous is None:
            bins = self.sav_binning()
            largest_bin = None
            max_bin_size = 0
            for bin_number, score_bin in enumerate(bins):
                if len(score_bin) > max_bin_size:
                    max_bin_size = len(score_bin)
                    largest_bin = bin_number

            median_neutral_bin = median(bins[largest_bin])
            median_synonymous = median_neutral_bin

        if self.mean_sav > median_synonymous: #larger score, stronger effect ...
            bottom_perc_median = median(sorted_scores[-perc_size:]) #...means, we need to take the right percentile of the sorted scores for a selection of scores with strong non-neutral effect
        else: #smaller score, stronger effect ...
            bottom_perc_median = median(sorted_scores[:perc_size]) #... means, we need to the left percentile

        if verbosity >= 1:
            print(f'Scaling SAV scores for {self.urn}')
            print(f'Min/Max score: {self.min_sav} / {self.max_sav}')
            print(f'Mean score: {self.mean_sav}')
            print(f'WT-like score: {median_synonymous}')
            print(f'Strong effect score: {bottom_perc_median}')

        if (bottom_perc_median - median_synonymous) == 0:
            return None

        for hgvs_pro in self.sav_scores:
            reported_score = self.sav_scores[hgvs_pro]
            scaled_score = score_scale_function(reported_score, median_synonymous, bottom_perc_median)
            self.scaled_sav_scores[hgvs_pro] = scaled_score

    def plot_sav_score_distribution(self, outfile):

        """
        Plots the histogram of the sav scores. Was majorly used for debugging purposes.
        Could be modified in the future to generate some additional output.

        Parameters
        ----------
        
        outfile
            Path to file, where the figure file should be written.
        """

        plt.rcParams.update({'figure.figsize':(7,5), 'figure.dpi':100})
        
        data = np.array(list(self.sav_scores.values()))

        plt.hist(data, bins=50)
        plt.gca().set(title=f'Score distribution for {self.urn}', ylabel='Frequency', xlabel = 'Effect score');
        plt.savefig(outfile)

@attr.s
class NewScoreSet(NewDataset):
    """
    This Class instantiates NewScoreSet and declares the variables present in NewScoreSet object.
    Attributes are set before posting a model instance
    """

    experiment: str = attr.ib(kw_only=True)

    # the following fields are optional
    meta_analysis_for: Optional[str] = attr.ib(kw_only=True, default=None)
    replaces: Optional[str] = attr.ib(kw_only=True, default=None)
    licence: Optional[str] = attr.ib(kw_only=True, default=None)
    data_usage_policy: Optional[str] = attr.ib(kw_only=True, default=None)

    # These can be strings or open filepaths
    score_data: Union[str, BinaryIO] = attr.ib(kw_only=True)
    count_data: Union[str, BinaryIO] = attr.ib(kw_only=True)
    meta_data: Union[str, BinaryIO] = attr.ib(kw_only=True)


@attr.s
class NewScoreSetRequest(APIObject):
    """
    This class instantiates NewScoreSetRequest and sets the fields for NewScoreSetRequest.
    Attributes are set before posting a model instance
    """

    scoreset: NewScoreSet = attr.ib(kw_only=True)
    target: NewTarget = attr.ib(kw_only=True)

    # These only need the ExternalIdentifier identifier field
    uniprot: SequenceOffset = attr.ib(kw_only=True, default=None)
    ensembl: SequenceOffset = attr.ib(kw_only=True, default=None)
    refseq: SequenceOffset = attr.ib(kw_only=True, default=None)

    reference_maps: List[ReferenceMap] = attr.ib(kw_only=True)

    def api_url() -> str:
        """
        Returns API endpoint
        """
        return "scoresets"

    def post_payload(self):
        """
        Use this to POST an instance of this class.
        data is converted to appropriate type and returned
        """
        json_dict, files = prepare_for_encoding(
            attr.asdict(
                self,
                filter=attrs_filter,
                retain_collection_types=True,
                value_serializer=attrs_serializer,
            )
        )
        return json_dict, files

