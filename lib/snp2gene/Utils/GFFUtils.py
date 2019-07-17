import os
import subprocess
import shutil
import csv
import pandas as pd
from pprint import pprint as pp


from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.GenomeSearchUtilClient import GenomeSearchUtil
from installed_clients.WorkspaceClient import Workspace


def tabix_index(filename, preset="gff", chrom="1", start="4", end="5", skip="0", comment="#"):
    """Call tabix to create an index for a bgzip-compressed file."""
    subprocess.Popen(['tabix', '-p', preset, '-s', chrom, '-b', start, '-e', end,
        '-S', skip, '-c', comment, filename])

def tabix_query(filename, chrom, start, end):
    """Call tabix and generate an array of strings for each line it returns."""
    query = f'{chrom}:{start}-{end}'
    process = subprocess.Popen(['tabix', '-f', filename, query], stdout=subprocess.PIPE)
    for line in process.stdout:
        yield line.decode('utf8').strip().split('\t')

def clean_tsv_data(data):
    data = str(data)
    data.replace('\n', '')
    return data

class GFFUtils:
    def __init__(self, config):
        self.callback_url = config['callback_url']
        self.shared_folder = config['scratch']
        self.ws_url = config['workspace-url']

        self.GFF_dir = os.path.join(self.shared_folder, 'GFF')

        if not os.path.isdir(self.GFF_dir):
            os.mkdir(self.GFF_dir)

        self.dfu = DataFileUtil(self.callback_url)
        self.gsu = GenomeSearchUtil(self.callback_url)
        self.wsc = Workspace(self.ws_url)

    def _prep_gff(self, gff_file):
        outfile = os.path.join(self.shared_folder, 'GFF', 'out.gff')
        sortcmd = f'(grep ^"#"  {gff_file}; grep -v ^"#" {gff_file} | sort -k1,1 -k4,4n)'
        with open(outfile, 'w') as o:
            p = subprocess.Popen(sortcmd, shell=True, stdout=o)
            p.wait()
            o.close()
        bgzip = subprocess.Popen(['bgzip', 'out.gff'], cwd=os.path.join(self.shared_folder, 'GFF'))
        bgzip.wait()
        outfile += '.gz'
        return outfile

    def _construct_gff_from_json(self, json, gff_file_path, contig_base_lengths):
        with open(gff_file_path, 'w') as f:
            for feature in json:
                if feature['feature_type'].strip().upper() == 'GENE':
                    end = int(feature['location'][0]['start'])+int(feature['location'][0]['length'])

                    metainfo = "ID="+feature['feature_id']

                    if feature['function']:
                        metainfo += ';FUNCTION='+feature['function']

                    contig_id = str(feature['location'][0]['contig_id'])
                    start = int(feature['location'][0]['start'])

                    # TODO: Fix Plink reassignment of Chr prefixes
                    try:
                        global_pos = int(contig_base_lengths[contig_id]) + start
                    except KeyError:
                        try:
                            global_pos = int(contig_base_lengths[contig_id.capitalize()]) + start
                        except KeyError:
                            try:
                                global_pos = int(contig_base_lengths['Chr'+str(contig_id)]) + start
                            except KeyError:
                                try:
                                    global_pos = int(contig_base_lengths['Chr0'+str(contig_id)]) + start
                                except KeyError:
                                    pp(contig_base_lengths)
                                    pp(contig_id)
                                    raise KeyError(e)

                    """
                    Remove ontology for now
                    if feature['ontology_terms']:
                        metainfo += ';ONTOLOGY('

                        for k, v in feature['ontology_terms'].items():
                            metainfo += str(k) + ',' + str(v) + ':'

                        metainfo = metainfo[:-1]  # remove trailing ;
                        metainfo += ')'
                    """

                    constructed_gff_line = str(feature['location'][0]['contig_id']) + '\t' + \
                                           'KBase\tgene\t' + \
                                           str(feature['location'][0]['start']) + '\t' + \
                                           str(end) + '\t.\t' + \
                                           str(feature['location'][0]['strand']) + '\t' + \
                                           str(global_pos) + '\t' + \
                                           str(metainfo) + '\n'
                    f.write(constructed_gff_line)
            f.close()
        if os.path.exists(gff_file_path):
            return gff_file_path
        else:
            raise FileNotFoundError('Unable to create GFF file form genome JSON.')

    def _process_tabix_results(self, queryresult):
        queryinfo = queryresult[8].split(';')
        if len(queryinfo) >= 2:
            extentsion = [clean_tsv_data(queryinfo[0][3:]), "NA", clean_tsv_data(queryinfo[1][9:])]
        elif len(queryinfo) is 1:
            extentsion = [clean_tsv_data(queryinfo[0][3:]), "NA", "NA"]
        else:
            extentsion = ['NA', 'NA', 'NA']
        return extentsion

    def find_geneid(self, row):
        tb = tabix_query(self.sorted_gff, row["CHR"], int(row["POS"]), int(row["POS"]))
        tbresult = next(tb, None)

        if tbresult is None:
            # do neighbor checking
            if int(row["POS"]) < 500:
                nstart = 0
            else:
                nstart = int(row["POS"]) - 500
            neigh_tb = tabix_query(self.sorted_gff, row["CHR"], nstart, int(row["POS"])+500)
            neigh_result = next(neigh_tb, None)

            if neigh_result is None:
                return pd.Series(['NA', 'NA', 'NA'], index=['GENEID','NEIGHBORGENE','FUNCTION'])
            else:
                nq = self._process_tabix_results(neigh_result)
                return pd.Series(nq, index=['GENEID', 'NEIGHBORGENE', 'FUNCTION'])
        else:
            q = self._process_tabix_results(tbresult)
            return pd.Series(q, index=['GENEID','NEIGHBORGENE','FUNCTION'])

    def annotate_GWAS_results(self, genome_ref, gwas_results_file):
        feature_num = self.gsu.search({'ref': genome_ref})['num_found']

        # get genome features for gff construction
        genome_features = self.gsu.search({
            'ref': genome_ref,
            'limit': feature_num,
            #'sort_by': [['feature_id', True]]
        })['features']

        assembly_ref = self.wsc.get_object_subset([{
            'included': ['/assembly_ref'],
            'ref': genome_ref
        }])[0]['data']['assembly_ref']

        # get assembly contigs for base length calculations
        assembly_contigs = self.wsc.get_object_subset([{
            'included': ['/contigs'],
            'ref': assembly_ref
        }])[0]['data']['contigs']

        contig_ids = list(assembly_contigs.keys())
        contig_ids.sort()

        contig_base_lengths = {}
        prev_length = 0

        for contig in contig_ids:
            contig_base_lengths[contig] = prev_length
            prev_length += assembly_contigs[contig]['length']

        gff_file = os.path.join(self.GFF_dir, 'constructed.gff')
        constructed_gff = self._construct_gff_from_json(genome_features, gff_file, contig_base_lengths)
        self.sorted_gff = self._prep_gff(constructed_gff)
        tabix_index(self.sorted_gff)

        gwas_results = pd.read_csv(gwas_results_file, sep='\t')

        gwas_results[['GENEID','NEIGHBORGENE','FUNCTION']] = \
            gwas_results.apply(self.find_geneid, axis=1)

        new_results_path = os.path.abspath(os.path.join(gwas_results_file, '..'))
        new_results_path = os.path.join(new_results_path, 'final_results.txt')

        gwas_results.to_csv(path_or_buf=new_results_path, sep='\t', index=False)

        return new_results_path
