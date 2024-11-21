import attr
import json
import os
from typing import Any, BinaryIO, Dict, List, Optional, Union
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

from .base import APIObject
from .dataset import Dataset
from .licence import Licence
from .target import Target
from .utils import attrs_filter, attrs_serializer, prepare_for_encoding, get_variant_type, median, score_scale_function, second_scale, disect_hgvs_single_pos
from fqfa.util.translate import translate_dna

from mavetools.models.sequence_retrieval import getUniprotSequence, get_refseq_sequences
from mavetools.models.utils import aac_to_hgvs_pro

@attr.s
class DatasetColumns:
    countColumns: List[str] = attr.ib(kw_only=True)
    scoreColumns: List[str] = attr.ib(kw_only=True)

@attr.s
class ScoreSet(APIObject, Dataset):
    """
    This class instantiates the ScoreSet object and declares the fields of that object.
    It inherits the attributes of APIObject and Dataset
    """

    urn: str = attr.ib(kw_only=True)
    license: Licence = attr.ib(kw_only=True)
    targetGenes: Target = attr.ib(kw_only=True)
    # optional attributes
    datasetColumns: Optional[DatasetColumns] = attr.ib(kw_only=True, default=None)
    replaces: Optional[str] = attr.ib(kw_only=True, default=None)
    #score_columns: List[str] = attr.ib(kw_only=True)
    #count_columns: List[str] = attr.ib(kw_only=True)
    # optional attributes
    previous_version: Optional[str] = attr.ib(kw_only=True, default=None)
    next_version: Optional[str] = attr.ib(kw_only=True, default=None)
    current_version: Optional[str] = attr.ib(kw_only=True, default=None)
    numVariants: int = attr.ib(kw_only=True, default=None)
    data_usage_policy: Optional[str] = attr.ib(kw_only=True, default=None)
    is_meta_analysis: Optional[bool] = attr.ib(kw_only=True, default=None)
    
    dataUsagePolicy: Optional[str] = attr.ib(kw_only=True, default=None)
    supersededScoreSetUrn: Optional[str] = attr.ib(kw_only=True, default=None)
    supersedingScoreSetUrn: Optional[str] = attr.ib(kw_only=True, default=None)
    metaAnalyzesScoreSetUrns: Optional[List[str]] = attr.ib(kw_only=True, default=None)
    metaAnalyzedByScoreSetUrns: Optional[List[str]] = attr.ib(kw_only=True, default=None)
    

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
        raw_seq = self.targetGenes[0]['targetSequence']['sequence']
        if self.targetGenes[0]['targetSequence']['sequenceType'] == 'protein':
            return raw_seq
        else:
            return translate_dna(raw_seq)[0]

    def get_full_sequence_info(self):
        if self.current_version == 'ProteinGym':
            full_seq = self.target['reference_sequence']['sequence']
            offset = 0
        else:
            id_dict = {}
            for identifier_struct in self.targetGenes[0]['externalIdentifiers']:
                dbname = identifier_struct['identifier']['dbName']
                offset = identifier_struct['offset']
                identifier = identifier_struct['identifier']['identifier']
                if dbname == 'UniProt':
                    id_dict['uniprot'] = offset, identifier
                elif dbname == 'Ensembl':
                    id_dict['ensembl'] = offset, identifier
                elif dbname == 'RefSeq':
                    id_dict['refseq'] = offset, identifier
            
            if 'uniprot' in id_dict:
                offset, identifier = id_dict['uniprot']
                full_seq = getUniprotSequence(identifier)
            #elif 'refseq' in id_dict:
            #    offset, identifier = id_dict['refseq']
            #    full_seq = get_refseq_sequences(identifier)
            #elif 'ensembl' in id_dict:
            #    print(f'No Uniprot, no RefSeq, but Ensembl ID: {self.urn}')
            #    full_seq = self.get_protein_sequence()
            #    offset = 0
            else:
                full_seq = self.get_protein_sequence()
                offset = 0
        return full_seq, offset

    def get_score_table_positions(self):
        exp_score_pos = None
        score_pos = None
        hgvs_pro_pos = 3
        hgvs_nt_pos = 1

        for pos, score_table_header in enumerate(self.datasetColumns['scoreColumns']):

            if score_table_header == 'score':
                score_pos = pos + 4
            if score_table_header == 'exp.score':
                exp_score_pos = pos + 4
        if exp_score_pos is not None:
            score_pos = exp_score_pos
        return hgvs_pro_pos, hgvs_nt_pos, score_pos

class ProteinGymScoreset(ScoreSet):
    def __init__(self, target_name, uniprot, seq, scoresetdata, offset = 0):
        self.target = {}
        self.target['uniprot'] = {}
        self.target['uniprot']['identifier'] = uniprot
        self.target['reference_sequence'] = {}
        self.target['reference_sequence']['sequence'] = seq
        self.target['reference_sequence']['sequence_type'] = 'protein'
        self.target['uniprot']['offset'] = offset
        self.target['name'] = target_name

        self.current_version = 'ProteinGym'
        self.scoresetdata = scoresetdata

def load_protein_gym_scores(dms_id, path_to_dms_file):

    f = open(path_to_dms_file)
    lines = f.readlines()
    f.close()

    score_dict = {}

    for line in lines[1:]:
        words = line.split(',')
        aac = words[0]
        if aac.count(':') > 0: #Consider only SAVs for the moment
            continue
        hgvs_pro = aac_to_hgvs_pro(aac)
        score = float(words[2])
        score_dict[hgvs_pro] = score

    ssd_object = ScoreSetData(dms_id, '', score_dict = score_dict)
    return ssd_object

class ScoreSetData:

    """
    This class represents the score tables for a corresponding scoresheet.
    Can also represent score tables for a corresponding experiment,
    in that case the scoretables from all scoresheets belonging to the experiment
    got aggregated (see aggregate_scoresetdata function of MLdataset class in ml_tools.py).
    """

    def __init__(self, urn, csv_formatted_data, hgvs_pro_pos = None, score_pos = None, hgvs_nt_pos = None, score_dict = None, verbosity = 1):

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
        if verbosity >= 2:
            print(f'Initialising ScoreSetData Object for {urn} {score_dict is None}')
        self.urn = urn
        if score_dict is None:
            self.score_dict = {}
        else:
            self.score_dict = score_dict

        self.nt_score_dict = {}
        self.sav_scores = {}
        self.synonymous_scores = {}
        self.nonsense_scores = {}
        self.mean_score = None
        self.min_score = None
        self.max_score = None
        self.mean_sav = None
        self.median_sav = None
        self.min_sav = None
        self.max_sav = None
        self.cardinality = 0
        self.sav_cardinality = 0
        self.sorted_sav_values = None
        self.mean_nonsense = None
        self.median_nonsense = None
        self.std_nonsense = None
        self.std_sav = None
        self.number_of_covered_positions = None
        self.normalizer = None

        #Parse the table

        lines = csv_formatted_data.split('\n')
        for line in lines:
            if line == '':
                continue
            if line[0] == '#':
                continue
            words = line.split(',')

            #In case their is a header row, parse the hgvs_pro, hgvs_nt, and effect column numbers
            if words[0] == 'accession':
                exp_score_pos = None
                for pos, word in enumerate(words):
                    if word == 'hgvs_pro':
                        hgvs_pro_pos = pos
                    if word == 'score':
                        score_pos = pos
                    if word == 'hgvs_nt':
                        hgvs_nt_pos = pos
                    if word == 'exp.score':
                        exp_score_pos = pos
                if exp_score_pos is not None:
                    score_pos = exp_score_pos
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

        #print(f'Parsed {len(self.score_dict)} number of scores for {self.urn}')

        #Get some basic information of the data


        self.pro_baseline = []
        self.ala_baseline = []

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
            elif variant_type == 'nonsense':
                self.nonsense_scores[hgvs_pro] = score
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
                if hgvs_pro[-3:] == 'Pro' and hgvs_pro[2:5] != 'Pro':
                    self.pro_baseline.append(score)

                if hgvs_pro[-3:] == 'Ala' and hgvs_pro[2:5] != 'Ala':
                    self.ala_baseline.append(score)

        #print(f'{len(self.sav_scores)} number of sav scores for {self.urn}')

        if len(self.pro_baseline) > 0:
            self.pro_baseline = sum(self.pro_baseline)/len(self.pro_baseline)
        else:
            self.pro_baseline = None

        if len(self.ala_baseline) > 0:
            self.ala_baseline = sum(self.ala_baseline)/len(self.ala_baseline)
        else:
            self.ala_baseline = None

        if self.cardinality > 0:
            self.mean_score = sum(self.score_dict.values()) / self.cardinality

        if self.sav_cardinality > 0:
            ssl = list(self.sav_scores.values())
            self.mean_sav = sum(ssl) / self.sav_cardinality
            self.std_sav = np.std(ssl)

        if len(self.nonsense_scores) > 0:
            nsl = list(self.nonsense_scores.values())
            self.median_nonsense = median(nsl)
            self.std_nonsense = np.std(nsl)


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


    def get_normalizer(self):
        if self.normalizer is not None:
            return self.normalizer
        sorted_scores = self.get_sorted_sav_values()
        bound = max([1,round(self.sav_cardinality/100)])
        self.normalizer = sorted_scores[self.sav_cardinality - bound] - sorted_scores[bound]
        return self.normalizer

    def get_number_of_covered_positions(self):
        if self.number_of_covered_positions is not None:
            return self.number_of_covered_positions
        covered_positions = set()
        for hgvs_pro in self.sav_scores:

            l, r, pos = disect_hgvs_single_pos(hgvs_pro[2:])
            covered_positions.add(pos)

        self.number_of_covered_positions = len(covered_positions)
        return self.number_of_covered_positions

    def sav_binning(self):

        """
        Creates a bin histogram of all sav scores.

        Returns
        -------

        bins
            The calculated histogram
        """

        sorted_scores = self.get_sorted_sav_values()
        n_of_bins = min([max([len(sorted_scores) // 100, 10]), len(sorted_scores)])
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


        bins = self.sav_binning()
        largest_bin = None
        highest_peak = None
        second_peak = None
        max_bin_size = 0
        highest_peak_size = 0
        second_peak_size = 0
        for bin_number, score_bin in enumerate(bins):
            if len(score_bin) > max_bin_size:
                max_bin_size = len(score_bin)
                largest_bin = bin_number

            #Is this a peak?
            peak = False
            if len(bins) == 1:
                peak = True
            elif bin_number == 0:
                if len(bins[1]) < len(score_bin):
                    peak = True
            elif bin_number == (len(bins) -1):
                if len(bins[bin_number - 1]) < len(score_bin):
                    peak = True
            else:
                if len(bins[bin_number - 1]) < len(score_bin) and len(bins[bin_number + 1]) < len(score_bin):
                    peak = True

            if peak:
                if len(score_bin) > highest_peak_size:
                    highest_peak_size = len(score_bin)
                    highest_peak = bin_number

                elif len(score_bin) > second_peak_size:
                    second_peak_size = len(score_bin)
                    second_peak = bin_number

        bi_modal = False
        if second_peak is not None:
            left_peak = min([highest_peak, second_peak])
            right_peak = max([highest_peak, second_peak])

            if second_peak_size < (0.4 * highest_peak_size):
                bi_modal = False
            else:
                for bin_between_peaks in bins[(left_peak+1):right_peak]:
                    if (2*len(bin_between_peaks)) < second_peak_size:
                        bi_modal = True

        if not bi_modal:
            median_neutral_bin = median(bins[largest_bin])
            binned_median_synonymous = median_neutral_bin

        else: #Currently we assume for a bimodal distribution that the neutral mutations got removed from the experiment (see urn:mavedb:00000071-a)
            median_highest_peak_bin = median(bins[highest_peak])
            median_second_peak_bin = median(bins[second_peak])
            binned_median_synonymous = 0.5 * (median_highest_peak_bin + median_second_peak_bin)

        if median_synonymous is None:
            calculated_by_binning = True
            median_synonymous = binned_median_synonymous
            if bi_modal:
                print(f'Bi-modal distribution for {self.urn}')
        else:
            calculated_by_binning = False
            bi_modal = None

        if self.mean_sav > median_synonymous: #larger score, stronger effect ...
            if self.pro_baseline is not None and self.ala_baseline is not None:
                if self.pro_baseline < self.ala_baseline:
                    print(f'Direction does not add up: {self.urn}')
                    bottom_perc_median = median(sorted_scores[:perc_size])
                    top_perc_median = median(sorted_scores[-perc_size:])
                else:
                    bottom_perc_median = median(sorted_scores[-perc_size:])
                    top_perc_median = median(sorted_scores[:perc_size])
            else:
                bottom_perc_median = median(sorted_scores[-perc_size:]) #...means, we need to take the right percentile of the sorted scores for a selection of scores with strong non-neutral effect
                top_perc_median = median(sorted_scores[:perc_size])
        else: #smaller score, stronger effect ...
            if self.pro_baseline is not None and self.ala_baseline is not None:
                if self.pro_baseline > self.ala_baseline:
                    print(f'Direction does not add up: {self.urn}')
                    bottom_perc_median = median(sorted_scores[-perc_size:])
                    top_perc_median = median(sorted_scores[:perc_size])
                else:
                    bottom_perc_median = median(sorted_scores[:perc_size])
                    top_perc_median = median(sorted_scores[-perc_size:])
            else:
                bottom_perc_median = median(sorted_scores[:perc_size]) #... means, we need to the left percentile
                top_perc_median = median(sorted_scores[-perc_size:])

        if verbosity >= 1:
            print(f'Scaling SAV scores for {self.urn}')
            print(f'Min/Max score: {self.min_sav} / {self.max_sav}')
            print(f'Mean score: {self.mean_sav}')
            print(f'WT-like score: {median_synonymous} (got calculated by binning: {calculated_by_binning}, binned median synonymous: {binned_median_synonymous}, bi_modal was {bi_modal})')
            print(f'Strong effect score: {bottom_perc_median}, Median nonsense: {self.median_nonsense}')
            print(f'Strong increasing effect score: {top_perc_median}')
            print(f'Pro baseline: {self.pro_baseline}')
            print(f'Ala baseline: {self.ala_baseline}')
            print(f'STD nonsense score: {self.std_nonsense}')
            print(f'STD missense score: {self.std_sav}')

        if (bottom_perc_median - median_synonymous) == 0:
            print(f'Bottom median equals to median synonymous: {self.urn}')
            return None

        for hgvs_pro in self.sav_scores:
            reported_score = self.sav_scores[hgvs_pro]
            scaled_score = score_scale_function(reported_score, median_synonymous, bottom_perc_median, top_perc_median)
            self.scaled_sav_scores[hgvs_pro] = scaled_score

        min_value = min(self.scaled_sav_scores.values())
        max_value = max(self.scaled_sav_scores.values())
        scaled_median = median(self.scaled_sav_scores.values())

        if scaled_median < 0.4:
            shift = 0.4 / scaled_median
        else:
            shift = None

        if shift is not None:
            print(f'Needed to perform a second scaling for {self.urn} with shift: {shift} (scaled median: {scaled_median})')
            for hgvs_pro in list(self.scaled_sav_scores.keys()):
                scaled_score = self.scaled_sav_scores[hgvs_pro]
                scaled_score = second_scale(scaled_score, min_value, max_value, shift)
                self.scaled_sav_scores[hgvs_pro] = scaled_score

        if verbosity >= 1:
            scaled_wt_score = score_scale_function(median_synonymous, median_synonymous, bottom_perc_median, top_perc_median)
            second_scaled_wt_score = second_scale(scaled_wt_score, min_value, max_value, shift)
            scaled_strong_score = score_scale_function(bottom_perc_median, median_synonymous, bottom_perc_median, top_perc_median)
            second_scaled_strong_score = second_scale(scaled_strong_score, min_value, max_value, shift)
            print(f'Scaled WT-like score: {scaled_wt_score}')
            print(f'Scaled strong effect score: {scaled_strong_score}')
            print(f'Scaled median effect score: {scaled_median}')
            print(f'Second scaled WT-like score: {second_scaled_wt_score}')
            print(f'Second scaled strong effect score: {second_scaled_strong_score}')

    def plot_sav_score_distribution(self, outfile, specific_scores = None):

        """
        Plots the histogram of the sav scores. Was majorly used for debugging purposes.
        Could be modified in the future to generate some additional output.

        Parameters
        ----------
        
        outfile
            Path to file, where the figure file should be written.

        specific_scores
            An optional score dictionary similar to the self.sav_scores object. Can be used for plotting other distributions, like the scaled ones.
        """

        plt.rcParams.update({'figure.figsize':(7,5), 'figure.dpi':100})
        
        if specific_scores is None:
            scores = self.sav_scores
        else:
            scores = specific_scores

        data = np.array(list(scores.values()))

        plt.hist(data, bins=50)
        plt.gca().set(title=f'Score distribution for {self.urn}', ylabel='Frequency', xlabel = 'Effect score');
        plt.savefig(outfile)


    def write_nonsense_tsv(self, outfile):
        lines = []
        for hgvs in self.nonsense_scores:
            score = self.nonsense_scores[hgvs]
            l, r, pos = disect_hgvs_single_pos(hgvs[2:])
            lines.append(f'{hgvs}\t{pos}\t{score}\n')

        f = open(outfile, 'w')
        f.write(''.join(lines))
        f.close()


    def calculate_nonsense_position_correlation(self):
        pos_score_tuples = []
        for hgvs in self.nonsense_scores:
            score = self.nonsense_scores[hgvs]
            l, r, pos = disect_hgvs_single_pos(hgvs[2:])
            pos_score_tuples.append((pos,score))
        pos_score_tuples = sorted(pos_score_tuples, key = lambda x: x[0])
        positions, scores = zip(*pos_score_tuples)
        pearson,pearson_p = stats.pearsonr(positions, scores)
        return pearson, pearson_p




