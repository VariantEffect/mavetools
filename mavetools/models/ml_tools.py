from mavetools.models.scoreset import ScoreSetData, ProteinGymScoreset, load_protein_gym_scores
from mavetools.models.utils import apply_offset_to_hgvs_pro, check_offset, median
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import metapub
import csv

class MlExperiment:

    """
    Class that represents one experiment.
    """

    def __init__(self, primary_id, scoreset_dict, scoreset, proteinGym_id = None, urn = None):

        """
        Initializes an instance of the ML experiment class.

        Parameters
        ----------

        urn
            MaveDB urn identifier of the experiment

        scoreset_dict
            A dictionary of all scoreset metdata objects belonging to the experiment.

        scoreset
            A representative metadata object for the experiment.
        """

        self.primary_id = primary_id
        self.urn = urn
        self.proteinGym_id = proteinGym_id
        self.scoreset_dict = scoreset_dict
        self.representative_scoreset_metadata = scoreset


    def get_name(self):
        if self.representative_scoreset_metadata.current_version == 'ProteinGym':
            return self.representative_scoreset_metadata.target['name']
        return self.representative_scoreset_metadata.targetGenes[0]['name'].replace(" ","_").replace('β','beta').replace(';','')

class MlDataset:

    """
    A class that represents a dataset of experiments determined for the usage in ML tools.
    """

    def __init__(self, experiments):

        """
        Initializes an instance of the ML dataset class.

        Parameters
        ----------

        experiments
            A dictionary mapping MaveDB urn identifiers to their corresponding ML experiment objects.
        """

        self.experiments = experiments

    def retrieve_data(self, client, verbosity = 0):

        """
        Uses a client object to download (or look up) the scoretables for all the scoresets in the dataset.

        Parameters
        ----------

        client
            A client object, see client.py for more information.
        """

        for experiment_urn in self.experiments:
            for scoreset_urn in self.experiments[experiment_urn].scoreset_dict:
                if verbosity >= 1:
                    print(f'Retrieving scoreset table {scoreset_urn} for experiment {experiment_urn}')
                csv_table = client.retrieve_score_table(scoreset_urn)
                hgvs_pro_pos, hgvs_nt_pos, score_pos = self.experiments[experiment_urn].scoreset_dict[scoreset_urn].get_score_table_positions()
                self.experiments[experiment_urn].scoreset_dict[scoreset_urn].scoresetdata = ScoreSetData(scoreset_urn, csv_table, hgvs_pro_pos = hgvs_pro_pos, score_pos = score_pos, hgvs_nt_pos = hgvs_nt_pos)


    def load_protein_gym(self, path_to_reference_file):

        proteingym_folder = path_to_reference_file.rsplit('/',1)[0]
        proteingym_folder = f'{proteingym_folder}/DMS_ProteinGym_substitutions'

        f = open(path_to_reference_file, 'r')
        lines = f.readlines()
        f.close()

        fetch = metapub.PubMedFetcher()
        doi_dict = {}
        for exp_urn in self.experiments:
            scoreset = self.experiments[exp_urn].representative_scoreset_metadata
            for publication_entry in scoreset.primaryPublicationIdentifiers:
                if publication_entry['publicationDoi'] is not None:
                    doi_dict[publication_entry['publicationDoi']] = exp_urn
                elif publication_entry['dbName'] == 'PubMed':
                    pubmed_id = publication_entry['identifier']
                    article = fetch.article_by_pmid(pubmed_id)
                    doi_dict[article.doi] = exp_urn

        for line in lines[1:]:
            words = [ '{}'.format(x) for x in list(csv.reader([line], delimiter=',', quotechar='"'))[0] ]

            dms_id = words[0]
            dms_filename = words[1]
            uniprot = words[2]
            wt_sequence = words[5]
            doi = words[16]
            if doi in doi_dict:
                print(f'Skipped ProteinGym entry: DOI already found in MaveDB: {dms_id} {doi_dict[doi]}')
                self.experiments[doi_dict[doi]].proteinGym_id = dms_id
                continue
            #else:
                #print(f'{dms_id} {doi}')
                
            path_to_dms_file = f'{proteingym_folder}/{dms_filename}'
            scoresetdata = load_protein_gym_scores(dms_id, path_to_dms_file)
            scoreset_obj = ProteinGymScoreset(uniprot, uniprot, wt_sequence, scoresetdata)
            self.experiments[dms_id] = MlExperiment(dms_id, {dms_id:scoreset_obj}, scoreset_obj, proteinGym_id = dms_id)


    def aggregate_scoresetdata(self, min_prot_size = None, min_coverage = None, min_len_coverage = None, reasonable_len_coverage = 100, std_filter = None, nonsense_std_filter = False, verbosity = 0):

        """
        For all experiments aggregates multiple scoresets corresponding to one experiment.
        Retrieves the full target protein sequence of the experiment.
        Checks all the offset values given in the metadata and maps all variations to the full target sequence.
        Effect values occupying the same position gets averaged.
        """

        del_list = []
        filtered_list = []

        for experiment_urn in self.experiments:
            if verbosity >= 1:
                print(f'Start aggregation of {experiment_urn}')
            experiment_scoresetdata = ScoreSetData(experiment_urn, '')
            experiment_scoresetdata.pro_baseline = []
            experiment_scoresetdata.ala_baseline = []

            broken = False
            reason = None

            for scoreset_urn in self.experiments[experiment_urn].scoreset_dict:
                if verbosity >= 1:
                    print(f'Processing data from scoreset {scoreset_urn} ...')
                scoreset = self.experiments[experiment_urn].scoreset_dict[scoreset_urn]
                scoresetdata = scoreset.scoresetdata

                seq, offset = scoreset.get_full_sequence_info()
                if seq is None:
                    reason = f'Filtered {scoreset_urn} due to no sequence available'
                    if verbosity >= 1:
                        print(reason)
                    broken = True
                    break

                if verbosity >= 2:
                    print(seq)

                if min_prot_size is not None:
                    if len(seq) < min_prot_size:
                        reason = f'Filtered {scoreset_urn} due to small prot size: {len(seq)} < {min_prot_size}'
                        if verbosity >= 1:
                            print(reason)
                        broken = True
                        break

                offset = check_offset(seq, offset, list(scoresetdata.sav_scores.keys()), scoreset_urn)

                if offset is None:
                    reason = f'Offset seems to be incorrect and could not be corrected for: {scoreset_urn}'
                    print(reason)
                    broken = True
                    break

                for hgvs_pro in scoresetdata.score_dict:
                    try:
                        corrected_hgvs_pro = apply_offset_to_hgvs_pro(hgvs_pro, offset)
                    except:
                        print(f'Couldnt apply offset to hgvs: {hgvs_pro} {offset}, {experiment_urn}')
                        corrected_hgvs_pro = apply_offset_to_hgvs_pro(hgvs_pro, offset)
                    if not corrected_hgvs_pro in experiment_scoresetdata.score_dict:
                        experiment_scoresetdata.score_dict[corrected_hgvs_pro] = []
                    experiment_scoresetdata.score_dict[corrected_hgvs_pro].append(scoresetdata.score_dict[hgvs_pro])

                for hgvs_pro in scoresetdata.sav_scores:
                    try:
                        corrected_hgvs_pro = apply_offset_to_hgvs_pro(hgvs_pro, offset)
                    except:
                        print(f'Couldnt apply offset to hgvs: {hgvs_pro} {offset}, {experiment_urn}')
                    if not corrected_hgvs_pro in experiment_scoresetdata.sav_scores:
                        experiment_scoresetdata.sav_scores[corrected_hgvs_pro] = []
                    score = scoresetdata.sav_scores[hgvs_pro]
                    experiment_scoresetdata.sav_scores[corrected_hgvs_pro].append(score)

                    if hgvs_pro[-3:] == 'Pro' and hgvs_pro[2:5] != 'Pro':
                        experiment_scoresetdata.pro_baseline.append(score)

                    if hgvs_pro[-3:] == 'Ala' and hgvs_pro[2:5] != 'Ala':
                        experiment_scoresetdata.ala_baseline.append(score)

                for hgvs_pro in scoresetdata.nonsense_scores:
                    try:
                        corrected_hgvs_pro = apply_offset_to_hgvs_pro(hgvs_pro, offset)
                    except:
                        print(f'Couldnt apply offset to hgvs: {hgvs_pro} {offset}, {experiment_urn}')
                    if not corrected_hgvs_pro in experiment_scoresetdata.nonsense_scores:
                        experiment_scoresetdata.nonsense_scores[corrected_hgvs_pro] = []
                    experiment_scoresetdata.nonsense_scores[corrected_hgvs_pro].append(scoresetdata.nonsense_scores[hgvs_pro])

                for hgvs_pro in scoresetdata.synonymous_scores:
                    try:
                        corrected_hgvs_pro = apply_offset_to_hgvs_pro(hgvs_pro, offset)
                    except:
                        corrected_hgvs_pro = hgvs_pro
                    if not corrected_hgvs_pro in experiment_scoresetdata.sav_scores:
                        experiment_scoresetdata.synonymous_scores[corrected_hgvs_pro] = []
                    experiment_scoresetdata.synonymous_scores[corrected_hgvs_pro].append(scoresetdata.synonymous_scores[hgvs_pro])

            if broken:
                del_list.append(experiment_urn)
                filtered_list.append((self.experiments[experiment_urn].get_name(), experiment_urn, reason))
                continue

            if len(experiment_scoresetdata.pro_baseline) > 0:
                experiment_scoresetdata.pro_baseline = sum(experiment_scoresetdata.pro_baseline)/len(experiment_scoresetdata.pro_baseline)
            else:
                experiment_scoresetdata.pro_baseline = None

            if len(experiment_scoresetdata.ala_baseline) > 0:
                experiment_scoresetdata.ala_baseline = sum(experiment_scoresetdata.ala_baseline)/len(experiment_scoresetdata.ala_baseline)
            else:
                experiment_scoresetdata.ala_baseline = None

            for hgvs_pro in experiment_scoresetdata.score_dict:
                experiment_scoresetdata.score_dict[hgvs_pro] = sum(experiment_scoresetdata.score_dict[hgvs_pro])/len(experiment_scoresetdata.score_dict[hgvs_pro])

            for hgvs_pro in experiment_scoresetdata.sav_scores:
                experiment_scoresetdata.sav_scores[hgvs_pro] = sum(experiment_scoresetdata.sav_scores[hgvs_pro])/len(experiment_scoresetdata.sav_scores[hgvs_pro])

            if len(experiment_scoresetdata.sav_scores) > 0:
                experiment_scoresetdata.std_sav = np.std(list(experiment_scoresetdata.sav_scores.values()))

            for hgvs_pro in experiment_scoresetdata.synonymous_scores:
                experiment_scoresetdata.synonymous_scores[hgvs_pro] = sum(experiment_scoresetdata.synonymous_scores[hgvs_pro])/len(experiment_scoresetdata.synonymous_scores[hgvs_pro])

            for hgvs_pro in experiment_scoresetdata.nonsense_scores:
                experiment_scoresetdata.nonsense_scores[hgvs_pro] = sum(experiment_scoresetdata.nonsense_scores[hgvs_pro])/len(experiment_scoresetdata.nonsense_scores[hgvs_pro])

            if len(experiment_scoresetdata.nonsense_scores) > 0:
                nsl = list(experiment_scoresetdata.nonsense_scores.values())
                experiment_scoresetdata.median_nonsense = median(nsl)
                experiment_scoresetdata.mean_nonsense = sum(nsl)/len(nsl)
                experiment_scoresetdata.std_nonsense = np.std(nsl)

            experiment_scoresetdata.mean_score = None
            if len(experiment_scoresetdata.score_dict) > 0:
                experiment_scoresetdata.min_score = min(experiment_scoresetdata.score_dict.values())
                experiment_scoresetdata.max_score = max(experiment_scoresetdata.score_dict.values())
            experiment_scoresetdata.mean_sav = None
            if len(experiment_scoresetdata.sav_scores) > 0:
                experiment_scoresetdata.min_sav = min(experiment_scoresetdata.sav_scores.values())
                experiment_scoresetdata.max_sav = max(experiment_scoresetdata.sav_scores.values())
            experiment_scoresetdata.cardinality = len(experiment_scoresetdata.score_dict)
            experiment_scoresetdata.sav_cardinality = len(experiment_scoresetdata.sav_scores)

            if min_coverage is not None:
                thresh = len(seq)*19*min_coverage
                if experiment_scoresetdata.sav_cardinality < thresh:
                    reason = f'Filtered {experiment_urn} due to variant coverage: {experiment_scoresetdata.sav_cardinality} ({experiment_scoresetdata.sav_cardinality/(len(seq)*19)}%) < {thresh} ({min_coverage}%)'
                    if verbosity >= 1:
                        print(reason)

                    filtered_list.append((self.experiments[experiment_urn].get_name(), experiment_urn, reason))
                    del_list.append(experiment_urn)
                    continue
            else:
                if verbosity >= 1:
                    print(f'{experiment_urn} variant coverage: {experiment_scoresetdata.sav_cardinality} ({experiment_scoresetdata.sav_cardinality/(len(seq)*19)}%)')

            number_of_covered_positions = experiment_scoresetdata.get_number_of_covered_positions()
            if min_len_coverage is not None:
                thresh = len(seq)*min_len_coverage
                if number_of_covered_positions < thresh and number_of_covered_positions < reasonable_len_coverage:
                    reason = f'Filtered {experiment_urn} due to length coverage: {number_of_covered_positions} ({number_of_covered_positions/(len(seq))}%) < {thresh} ({min_len_coverage}%)'
                    if verbosity >= 1:
                        print(reason)
                    del_list.append(experiment_urn)
                    filtered_list.append((self.experiments[experiment_urn].get_name(), experiment_urn, reason))
                    continue
            else:
                if verbosity >= 1:
                    print(f'{experiment_urn} length coverage: {number_of_covered_positions} ({number_of_covered_positions/(len(seq))}%)')

            if nonsense_std_filter:
                if experiment_scoresetdata.std_nonsense is not None:
                    if verbosity >= 1:
                        if len(experiment_scoresetdata.nonsense_scores) == 1:
                            print("Can't calculate Pearson's correlation coefficient between nonsense mutation effect scores and corresponding position number, due there is only one value")
                        else:
                            pearson_corr, pearson_corr_p_value = experiment_scoresetdata.calculate_nonsense_position_correlation()
                            print(f"{experiment_urn}: Pearson's correlation coefficient between nonsense mutation effect scores and corresponding position number: {pearson_corr} (p-value: {pearson_corr_p_value})")
                    if experiment_scoresetdata.std_sav is not None:
                        if experiment_scoresetdata.std_nonsense*5 > experiment_scoresetdata.std_sav:
                            reason = f'Filtered {experiment_urn} due to high nonsense deviation: {experiment_scoresetdata.std_nonsense} STD nonsense score versus {experiment_scoresetdata.std_sav} STD missense score'
                            if verbosity >= 1:
                                print(reason)
                            del_list.append(experiment_urn)
                            filtered_list.append((self.experiments[experiment_urn].get_name(), experiment_urn, reason))
                            continue

            if std_filter is not None:
                normalizer = experiment_scoresetdata.get_normalizer()
                std_missense = experiment_scoresetdata.std_sav
                normalized_std_missense = std_missense/normalizer

                if normalized_std_missense < std_filter:
                    reason = f'Filtered {experiment_urn} due to low normalized missense STD: {normalized_std_missense} normalized missense STD versus < {std_filter}'
                    if verbosity >= 1:
                        print(reason)
                    del_list.append(experiment_urn)
                    filtered_list.append((self.experiments[experiment_urn].get_name(), experiment_urn, reason))
                    continue

            if experiment_scoresetdata.cardinality > 0:
                experiment_scoresetdata.mean_score = sum(experiment_scoresetdata.score_dict.values()) / experiment_scoresetdata.cardinality

            if experiment_scoresetdata.sav_cardinality > 0:
                experiment_scoresetdata.mean_sav = sum(experiment_scoresetdata.sav_scores.values()) / experiment_scoresetdata.sav_cardinality
                experiment_scoresetdata.median_sav = median(experiment_scoresetdata.sav_scores.values())

            self.experiments[experiment_urn].experiment_scoresetdata = experiment_scoresetdata

        for experiment_urn in del_list:
            del self.experiments[experiment_urn]

        return filtered_list

    def write_filtered_entries(self, filtered_list, outfile):

        outlines = []

        for name, experiment_urn, reason in filtered_list:

            outlines.append(f'{experiment_urn}\t{name}\t{reason}\n')

        f = open(outfile,'w')
        f.write(''.join(outlines))
        f.close()


    def scale_all_savs(self, verbosity = 0):

        """
        Applies to effect value scaling to all SAV effect scores in the experiment.

        Parameters
        ----------

        verbosity
            Verbosity level of the function.
        """

        for experiment_urn in self.experiments:
            if len(self.experiments[experiment_urn].experiment_scoresetdata.sav_scores) == 0:
                continue
            self.experiments[experiment_urn].experiment_scoresetdata.scale_sav_data(verbosity = verbosity)

    def plot_sav_score_distribution(self, outfile):
        all_scores = []
        for experiment_urn in self.experiments:
            all_scores += self.experiments[experiment_urn].experiment_scoresetdata.scaled_sav_scores.values()

        plt.rcParams.update({'figure.figsize':(7,5), 'figure.dpi':100})

        data = np.array(all_scores)

        plt.hist(data, bins=50)
        plt.gca().set(title=f'Score distribution', ylabel='Frequency', xlabel = 'Effect score');
        plt.savefig(outfile)


    def write_dataset_statistics(self, statfile):
        headers = [
            'Experiment URN',
            'Protein Length',
            'Dataset Cardinality',
            'Coverage',
            'Length Coverage',
            'Value Range',
            'Normalizer',
            'Missense Mean',
            'Missense Median',
            'Nonsense Mean',
            'Nonsense Median',
            'Missense STD',
            'Normalized Missense STD',
            'Length Normalized Missense STD'
        ]
        header = "\t".join(headers)
        outlines = [f'{header}\n']

        for experiment_urn in self.experiments:

            name = self.experiments[experiment_urn].representative_scoreset_metadata.target['name'].replace(" ","_").replace('β','beta')

            seq, offset = self.experiments[experiment_urn].representative_scoreset_metadata.get_full_sequence_info()
            prot_len = len(seq)
            data_cardinality = self.experiments[experiment_urn].experiment_scoresetdata.sav_cardinality
            coverage = data_cardinality/(prot_len*19)
            number_of_covered_positions = self.experiments[experiment_urn].experiment_scoresetdata.get_number_of_covered_positions()
            length_coverage = number_of_covered_positions/prot_len

            min_val = self.experiments[experiment_urn].experiment_scoresetdata.min_sav
            max_val = self.experiments[experiment_urn].experiment_scoresetdata.max_sav

            if min_val == max_val:
                print(f'Min val equals to max val for {experiment_urn}')
                continue

            normalizer = self.experiments[experiment_urn].experiment_scoresetdata.get_normalizer()

            std_missense = self.experiments[experiment_urn].experiment_scoresetdata.std_sav

            normalized_std_missense = std_missense/normalizer

            length_normalized_std_missense = std_missense/(max_val-min_val)

            val_range = f'[{min_val} : {max_val}]'

            missense_mean = self.experiments[experiment_urn].experiment_scoresetdata.mean_sav
            missense_median = self.experiments[experiment_urn].experiment_scoresetdata.median_sav
            nonsense_mean = self.experiments[experiment_urn].experiment_scoresetdata.mean_nonsense
            nonsense_median = self.experiments[experiment_urn].experiment_scoresetdata.median_nonsense

            line = '\t'.join([
                f'{name}_{experiment_urn}',
                f'{prot_len}',
                f'{data_cardinality}',
                f'{coverage}',
                f'{length_coverage}',
                f'{val_range}',
                f'{normalizer}',
                f'{missense_mean}',
                f'{missense_median}',
                f'{nonsense_mean}',
                f'{nonsense_median}',
                f'{std_missense}',
                f'{normalized_std_missense}',
                f'{length_normalized_std_missense}'
            ])

            outlines.append(f'{line}\n')

        f = open(statfile, 'w')
        f.write(''.join(outlines))
        f.close()


    def write_scaled_sav_fasta(self, outfile, sequences_only_file = None, write_unscaled = False):

        """
        Writes scaled SAV effect scores together with their full target protein sequences to a fasta file.
        The fasta file has a non-standard format:

            >[Sequence identifier]
            <[hgvs_pro identifier] #scaled_effect:[scaled effect value]
            <[hgvs_pro identifier] #scaled_effect:[scaled effect value]
            <[hgvs_pro identifier] #scaled_effect:[scaled effect value]
            ...
            [The sequence]

        Parameters
        ----------

        outfile
            Path to where the fasta file should be written.

        sequences_only_file
            When not None, gives a path where a classical fasta gets written omitting the variation information.
        """

        fasta_lines = []
        seq_only_lines = []
        for experiment_urn in self.experiments:
            seq, offset = self.experiments[experiment_urn].representative_scoreset_metadata.get_full_sequence_info()

            name = self.experiments[experiment_urn].get_name()

            try:
                scaled_sav_scores = self.experiments[experiment_urn].experiment_scoresetdata.scaled_sav_scores
            except:
                print(f'Couldnt find scaled sav scores for {experiment_urn}')
                continue
            if write_unscaled:
                scaled_sav_scores = self.experiments[experiment_urn].experiment_scoresetdata.sav_scores
            if len(scaled_sav_scores) == 0:
                print(f'Scaled sav scores empty for {experiment_urn}')
                continue

            headline = f'>{name}_{experiment_urn}\t<{self.experiments[experiment_urn].urn},{self.experiments[experiment_urn].proteinGym_id}\n'
            fasta_lines.append(headline)
            if sequences_only_file is not None:
                seq_only_lines.append(headline)

            for hgvs_pro in scaled_sav_scores:

                variant_line = f'<{hgvs_pro} #scaled_effect:{scaled_sav_scores[hgvs_pro]}\n'
                fasta_lines.append(variant_line)

            fasta_lines.append(f'{seq}\n')
            if sequences_only_file is not None:
                seq_only_lines.append(f'{seq}\n')

        f = open(outfile, 'w')
        f.write(''.join(fasta_lines))
        f.close()

        if sequences_only_file is not None:
            f = open(sequences_only_file, 'w')
            f.write(''.join(seq_only_lines))
            f.close()

