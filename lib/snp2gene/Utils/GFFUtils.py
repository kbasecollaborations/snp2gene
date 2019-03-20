import os
import subprocess
import shutil
import csv
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
        self.GFF_dir = os.path.join(self.shared_folder, 'GFF')

        if not os.path.isdir(self.GFF_dir):
            os.mkdir(self.GFF_dir)

        self.dfu = DataFileUtil(self.callback_url)
        self.gsu = GenomeSearchUtil(self.callback_url)
        # TODO: replace hard coded url with config['workspace_url']
        self.wsc = Workspace("https://appdev.kbase.us/services/ws")

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

                    global_pos = int(contig_base_lengths[contig_id]) + start

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
        queryinfo = queryresult[7].split(';')
        if len(queryinfo) >= 2:
            extentsion = [clean_tsv_data(queryinfo[0][3:]), "NA", clean_tsv_data(queryinfo[1][9:])]
        elif len(queryinfo) is 1:
            extentsion = [clean_tsv_data(queryinfo[0][3:]), "NA", "NA"]
        else:
            extentsion = ['NA', 'NA', 'NA']
        return extentsion

    def annotate_GWAS_results(self, genome_ref, gwas_results_file):
        feature_num = self.gsu.search({'ref': genome_ref})['num_found']

        # get genome features for gff construction
        genome_features = self.gsu.search({
            'ref': genome_ref,
            'limit': feature_num,
            'sort_by': [['feature_id', True]]
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

        contig_base_lengths = {}
        prev_length = 0

        for contig in assembly_contigs:
            contig_base_lengths[contig] = prev_length
            prev_length += assembly_contigs[contig]['length']

        gff_file = os.path.join(self.GFF_dir, 'constructed.gff')
        constructed_gff = self._construct_gff_from_json(genome_features, gff_file, contig_base_lengths)
        sorted_gff = self._prep_gff(constructed_gff)
        tabix_index(sorted_gff)

        new_results_file = []

        with open(gwas_results_file, 'r') as gwasresults:
            gwasreader = csv.reader(gwasresults, delimiter='\t')
            next(gwasreader)  # skip headers
            for result in gwasreader:
                tb = tabix_query(sorted_gff, 'Chr'+result[1], int(result[2]), int(result[2]))
                tbquery = next(tb, None) # if there is no first object in generator, set to None
                if tbquery is not None:
                    query_result = self._process_tabix_results(tbquery)
                    result.extend(query_result)
                else:
                    tb_neighbors = tabix_query(sorted_gff, 'Chr' + result[1], int(result[2])-1000,
                                                int(result[2])+1000)
                    tbquery_neighbors = next(tb_neighbors, None)
                    if tbquery_neighbors is not None:
                        query_neighbor_result = self._process_tabix_results(tbquery_neighbors)
                        result.extend(query_neighbor_result)
                    else:
                        result.extend(['NA', 'NA', 'NA'])

                new_results_file.append(result)

            gwasresults.close()

        new_results_headers = "SNP\tCHR\tBP\tP\tPOS\tGENEID\tNEIGHBORGENE\tFUNCTION\n"

        if not os.access(gwas_results_file, os.W_OK):
            shutil.copyfile(gwas_results_file, os.path.join(self.shared_folder, os.path.basename(gwas_results_file)))
            gwas_results_file = os.path.join(self.shared_folder, os.path.basename(gwas_results_file))

        with open(gwas_results_file, 'w') as r:
            r.write(new_results_headers)
            for line in new_results_file:
                r.write('\t'.join(line)+'\n')
            r.close()

        return gwas_results_file
