from mavetools.models.scoreset import ScoreSetData
from mavetools.models.utils import apply_offset_to_hgvs_pro, check_offset

class MlExperiment:
    def __init__(self, urn, scoreset_dict, scoreset):
        self.urn = urn
        self.scoreset_dict = scoreset_dict
        self.representative_scoreset_metadata = scoreset


class MlDataset:
    def __init__(self, experiments):
        self.experiments = experiments

    def retrieve_data(self, client):
        for experiment_urn in self.experiments:
            for scoreset_urn in self.experiments[experiment_urn].scoreset_dict:
                csv_table = client.retrieve_score_table(scoreset_urn)
                hgvs_pro_pos, hgvs_nt_pos, score_pos = self.experiments[experiment_urn].scoreset_dict[scoreset_urn].get_score_table_positions()
                self.experiments[experiment_urn].scoreset_dict[scoreset_urn].scoresetdata = ScoreSetData(scoreset_urn, csv_table, hgvs_pro_pos = hgvs_pro_pos, score_pos = score_pos, hgvs_nt_pos = hgvs_nt_pos)

    def aggregate_scoresetdata(self):
        for experiment_urn in self.experiments:
            
            experiment_scoresetdata = ScoreSetData(experiment_urn, '')

            for scoreset_urn in self.experiments[experiment_urn].scoreset_dict:
                scoreset = self.experiments[experiment_urn].scoreset_dict[scoreset_urn]
                scoresetdata = scoreset.scoresetdata

                seq, offset = scoreset.get_full_sequence_info()

                offset = check_offset(seq, offset, list(scoresetdata.sav_scores.keys()), scoreset_urn)

                if offset is None:
                    print(f'Offset seems to be incorrect and could not be corrected for: {scoreset_urn}')
                    break

                for hgvs_pro in scoresetdata.score_dict:
                    if not hgvs_pro in experiment_scoresetdata.score_dict:
                        experiment_scoresetdata.score_dict[hgvs_pro] = []
                    experiment_scoresetdata.score_dict[hgvs_pro].append(scoresetdata.score_dict[hgvs_pro])

                for hgvs_pro in scoresetdata.sav_scores:
                    try:
                        corrected_hgvs_pro = apply_offset_to_hgvs_pro(hgvs_pro, offset)
                    except:
                        print(f'Couldnt apply offset to hgvs: {hgvs_pro} {offset}, {experiment_urn}')
                    if not corrected_hgvs_pro in experiment_scoresetdata.sav_scores:
                        experiment_scoresetdata.sav_scores[corrected_hgvs_pro] = []
                    experiment_scoresetdata.sav_scores[corrected_hgvs_pro].append(scoresetdata.sav_scores[hgvs_pro])

                for hgvs_pro in scoresetdata.synonymous_scores:
                    if not hgvs_pro in experiment_scoresetdata.sav_scores:
                        experiment_scoresetdata.synonymous_scores[hgvs_pro] = []
                    experiment_scoresetdata.synonymous_scores[hgvs_pro].append(scoresetdata.synonymous_scores[hgvs_pro])

            for hgvs_pro in experiment_scoresetdata.score_dict:
                experiment_scoresetdata.score_dict[hgvs_pro] = sum(experiment_scoresetdata.score_dict[hgvs_pro])/len(experiment_scoresetdata.score_dict[hgvs_pro])

            for hgvs_pro in experiment_scoresetdata.sav_scores:
                experiment_scoresetdata.sav_scores[hgvs_pro] = sum(experiment_scoresetdata.sav_scores[hgvs_pro])/len(experiment_scoresetdata.sav_scores[hgvs_pro])

            for hgvs_pro in experiment_scoresetdata.synonymous_scores:
                experiment_scoresetdata.synonymous_scores[hgvs_pro] = sum(experiment_scoresetdata.synonymous_scores[hgvs_pro])/len(experiment_scoresetdata.synonymous_scores[hgvs_pro])

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

            if experiment_scoresetdata.cardinality > 0:
                experiment_scoresetdata.mean_score = sum(experiment_scoresetdata.score_dict.values()) / experiment_scoresetdata.cardinality

            if experiment_scoresetdata.sav_cardinality > 0:
                experiment_scoresetdata.mean_sav = sum(experiment_scoresetdata.sav_scores.values()) / experiment_scoresetdata.sav_cardinality

            self.experiments[experiment_urn].experiment_scoresetdata = experiment_scoresetdata

    def scale_all_savs(self, verbosity = 0):
        for experiment_urn in self.experiments:
            if len(self.experiments[experiment_urn].experiment_scoresetdata.sav_scores) == 0:
                continue
            self.experiments[experiment_urn].experiment_scoresetdata.scale_sav_data(verbosity = verbosity)

    def write_scaled_sav_fasta(self, outfile, sequences_only_file = None):
        fasta_lines = []
        seq_only_lines = []
        for experiment_urn in self.experiments:
            seq, offset = self.experiments[experiment_urn].representative_scoreset_metadata.get_full_sequence_info()

            name = self.experiments[experiment_urn].representative_scoreset_metadata.target['name'].replace(" ","_").replace('β','beta')

            try:
                scaled_sav_scores = self.experiments[experiment_urn].experiment_scoresetdata.scaled_sav_scores
            except:
                continue
            if len(scaled_sav_scores) == 0:
                continue

            headline = f'>{name}_{experiment_urn}\n'
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

def create_sav_ml_dataset(client, scoreset_dict, outfile):
    fasta_lines = []
    for urn in scoreset_dict:
        csv_table = client.retrieve_score_table(urn)
        scoresetdata = ScoreSetData(csv_table)
        normalized_sav_scores = scoresetdata.normalize_sav_data()
        seq = scoreset_dict[urn].get_protein_sequence()

        headline = f'>{scoreset_dict[urn].target.name}\n'
        fasta_lines.append(headline)
        for hgvs_pro in normalized_sav_scores:
            variant_line = f'<{hgvs_pro} {normalized_sav_scores[hgvs_pro]}\n'
            fasta_lines.append(variant_line)

        fasta_lines.append(f'{seq}\n')

    f = open(outfile, 'w')
    f.write(''.join(fasta_lines))
    f.close()
