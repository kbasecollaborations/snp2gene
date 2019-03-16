import os
import subprocess
import json
import csv
from pprint import pprint as pp
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.GenomeFileUtilClient import GenomeFileUtil
from installed_clients.GenomeSearchUtilClient import GenomeSearchUtil


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
        self.gfu = GenomeFileUtil(self.callback_url)
        self.gsu = GenomeSearchUtil(self.callback_url)

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

    def _construct_gff_from_json(self, json, gff_file_path):
        with open(gff_file_path, 'w') as f:
            for feature in json:
                if feature['feature_type'].strip().upper() == 'GENE':
                    end = int(feature['location'][0]['start'])+int(feature['location'][0]['length'])

                    metainfo = "ID="+feature['feature_id']

                    if feature['function']:
                        metainfo += ';FUNCTION='+feature['function']

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
                                           metainfo + '\n'
                    f.write(constructed_gff_line)
            f.close()

        if os.path.exists(gff_file_path):
            return gff_file_path
        else:
            raise FileNotFoundError('Unable to create GFF file form genome JSON.')

    def annotate_GWAS_results(self, genome_ref, association_ref):
        """
        assoc_results = self.dfu.get_objects({'object_refs': [association_ref]})['data'][0]['data'][
            'association_details']

            
        feature_num = self.gsu.search({'ref': genome_ref})['num_found']
        genome_features = self.gsu.search({
            'ref': genome_ref,
            'limit': feature_num,
            'sort_by': [['feature_id', True]]
        })['features']
        """
        with open('/kb/module/work/features.json', 'r') as h:
            genome_features = json.load(h)

        gff_file = os.path.join(self.GFF_dir, 'constructed.gff')
        constructed_gff = self._construct_gff_from_json(genome_features, gff_file)
        sorted_gff = self._prep_gff(constructed_gff)
        print(sorted_gff)
        tabix_index(sorted_gff)

        gwas_results_file = '/kb/module/work/snpdataFLC.tsv'
        new_results_file = []

        # first pass fill in snps that reside in genes
        with open(gwas_results_file, 'r') as gwasresults:
            gwasreader = csv.reader(gwasresults, delimiter='\t')
            next(gwasreader)  # skip headers
            for result in gwasreader:
                tb = tabix_query(sorted_gff, 'Chr'+result[1], int(result[2]), int(result[2]))
                tbquery = next(tb, None) # if there is no first object in generator, set to None
                if tbquery is not None:
                    queryinfo = tbquery[7].split(';')
                    if len(queryinfo) >= 2:
                        result.extend([clean_tsv_data(queryinfo[0][3:]), "NA",
                                       clean_tsv_data(queryinfo[1][9:])])
                    elif len(queryinfo) == 1:
                        result.extend([clean_tsv_data(queryinfo[0][3:]), "NA", "NA"])
                    else:
                        # csv addition columns are GENEID NEIGHBORGENE FUNCTION
                        result.extend(['NA', 'NA', 'NA'])
                else:
                    tb_neighbors = tabix_query(sorted_gff, 'Chr' + result[1], int(result[2])-1000,
                                                int(result[2])+1000)
                    tbquery_neighbors = next(tb_neighbors, None)

                    if tbquery_neighbors is not None:
                        queryinfo_neighbor = tbquery_neighbors[7].split(';')
                        if len(queryinfo_neighbor) >= 2:
                            result.extend([clean_tsv_data(queryinfo_neighbor[0][3:]), "NA",
                                           clean_tsv_data(queryinfo_neighbor[1][9:])])
                        elif len(queryinfo_neighbor) == 1:
                            result.extend([clean_tsv_data(queryinfo_neighbor[0][3:]), "NA", "NA"])
                        else:
                            # csv addition columns are GENEID NEIGHBORGENE FUNCTION
                            result.extend(['NA', 'NA', 'NA'])
                    else:
                        result.extend(['NA', 'NA', 'NA'])

                new_results_file.append(result)

        new_results_file_path = os.path.join(self.shared_folder, 'newSNPdata.tsv')
        new_results_headers = "SNP\tCHR\tBP\tP\tPOS\tGENEID\tNEIGHBORGENE\tFUNCTION\n"

        with open(new_results_file_path, 'w') as r:
            r.write(new_results_headers)
            for line in new_results_file:
                r.write('\t'.join(line)+'\n')
            r.close()

        return new_results_file
