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

class GFFUtils2:
    def __init__(self, config):
        self.callback_url = config['callback_url']
        self.shared_folder = config['scratch']
        #self.shared_folder = "/kb/module/work"
        self.ws_url = config['workspace-url']

        self.dfu = DataFileUtil(self.callback_url)
        self.gsu = GenomeSearchUtil(self.callback_url)
        self.wsc = Workspace(self.ws_url)

    def _prep_gff(self, gff_file):
        outfile = os.path.join(self.genome_dir, 'out.gff')
        sortcmd = f'(grep ^"#"  {gff_file}; grep -v ^"#" {gff_file} | sort -k1,1 -k4,4n)'

        with open(outfile, 'w') as o:
            p = subprocess.Popen(sortcmd, shell=True, stdout=o)
            out, err = p.communicate()
            o.close()

        bgzip = subprocess.Popen(['bgzip', 'out.gff'], cwd=self.genome_dir)
        out2, err2 = bgzip.communicate()

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
            extension = [clean_tsv_data(queryinfo[0][3:]), "NA", clean_tsv_data(queryinfo[1][9:])]
        elif len(queryinfo) is 1:
            extension = [clean_tsv_data(queryinfo[0][3:]), "NA", "NA"]
        else:
            extension = ['NA', 'NA', 'NA']
        return extension

    def find_gene_info(self, row):
        tb = tabix_query(self.sorted_gff, row["CHR"], int(row["POS"]), int(row["POS"]))
        tbresult = next(tb, None)
        if tbresult is None:
            tb2 = tabix_query(self.sorted_gff, 'chr' + row["CHR"], int(row["POS"]), int(row["POS"]))
            tbresult2 = next(tb2, None)
            if tbresult2 is None:
                tb3 = tabix_query(self.sorted_gff, 'chr0' + row["CHR"], int(row["POS"]), int(row["POS"]))
                tbresult3 = next(tb3, None)
                if tbresult3 is None:
                    if int(row["POS"]) < 500:
                        nstart = 0
                    else:
                        nstart = int(row["POS"]) - 500

                    neigh_tb = tabix_query(self.sorted_gff, row["CHR"], nstart, int(row["POS"]) + 500)
                    neigh_result = next(neigh_tb, None)

                    if neigh_result is None:
                        return pd.Series(['NA', 'NA', 'NA'], index=['GENEID', 'NEIGHBORGENE', 'FUNCTION'])
                    else:
                        nq = self._process_tabix_results(neigh_result)
                        return pd.Series([nq[1], nq[0], nq[2]], index=['GENEID', 'NEIGHBORGENE', 'FUNCTION'])
                else:
                    q3 = self._process_tabix_results(tbresult3)
                    return pd.Series(q3, index=['GENEID', 'NEIGHBORGENE', 'FUNCTION'])
            else:
                q2 = self._process_tabix_results(tbresult2)
                return pd.Series(q2, index=['GENEID', 'NEIGHBORGENE', 'FUNCTION'])
        else:
            q = self._process_tabix_results(tbresult)
            return pd.Series(q, index=['GENEID', 'NEIGHBORGENE', 'FUNCTION'])

    def get_gwas_result_file(self, association_ref, association_name, p_value):
        #association_obj = self.dfu.get_objects({'object_refs': [association_ref]})['data'][0]['data']['data']
        association_obj = self.dfu.get_objects({'object_refs': [association_ref]})['data'][0]
        association_results = association_obj['data']["association_details"][0]["association_results"]
        result = "CHR\tSNP\tPOS\tP\tBP\n"
        for variation in association_results:
            if (float(variation[3]) > float(p_value)):
                continue
            result += str(variation[0]) + "\t" 
            result +=  str(variation[1]) + "\t" 
            result +=  str(variation[2]) + "\t" 
            result +=   str(variation[3]) + "\t"
            result +=   str(variation[2]) + "\n"
        filepath = os.path.join(self.genome_dir, association_name)
        with open(filepath, "w") as file1: 
            file1.write(result) 
        return (filepath)

    def build_featureset(self, filepath, genome_ref, description, workspace_name, association_name, prefix):
      gene_ids = dict()
      element_ordering = list()
      elements = dict()
      skip_words = ["GENEID", "NEIGHBORGENE", "NA"]
      with open(filepath, 'r') as reader:
          for line in reader:
              fields = line.split("\t")
              condition1 = fields[5] not in skip_words
              condition2 = fields[5] not in elements
              condition3 = fields[6] not in skip_words
              condition4 = fields[6] not in elements
              if condition1 and condition2:
                  element_ordering.append(fields[5])
                  elements[fields[5]] = [genome_ref]
              if condition3 and condition4:
                  element_ordering.append(fields[6])
                  elements[fields[6]] = [genome_ref]
      featureset = dict()
      featureset['description'] = description
      featureset['element_ordering'] = element_ordering
      featureset['elements'] = elements
      ws_id = self.dfu.ws_name_to_id(workspace_name)
      featureset_obj_name = prefix + str(association_name)

      save_info = self.dfu.save_objects( { 'id': ws_id, 
                                            'objects': [ {'type': 'KBaseCollections.FeatureSet', 
                                                          'data': featureset, 
                                                          'name': featureset_obj_name}]})[0]
      obj_ref  = "{0}/{1}/{2}".format( save_info[6], save_info[0], save_info[4] )   
      return obj_ref         


   
    def annotate_GWAS_results(self, genome_ref, association_ref, workspace_name, prefix, p_value):
         
        #TODO: Send outfile to prep gff function inseted of hardcord
        #TODO: Removed hard coded stuff and create new directory for each test function
        self.genome_dir_name = "_".join(genome_ref.split("/"))
        self.genome_dir = os.path.join(self.shared_folder, self.genome_dir_name)
        if not os.path.isdir(self.genome_dir):
            os.mkdir(self.genome_dir)
        sorted_gff_path = os.path.join(self.genome_dir, 'out.gff.gz')
        self.sorted_gff = sorted_gff_path

        if  not os.path.exists(sorted_gff_path):
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

            gff_file = os.path.join(self.genome_dir, 'constructed.gff')
            constructed_gff = self._construct_gff_from_json(genome_features, gff_file, contig_base_lengths)
            self.sorted_gff = self._prep_gff(constructed_gff)
            tabix_index(self.sorted_gff)

        obj_info = self.wsc.get_object_info3({"objects": [{"ref": association_ref}]})
        association_name =obj_info["infos"][0][1]


        gwas_results_file = self.get_gwas_result_file(association_ref, association_name, p_value)

        gwas_results = pd.read_csv(gwas_results_file, sep='\t')

        gwas_results[['GENEID', 'NEIGHBORGENE', 'FUNCTION']] = \
           gwas_results.apply(self.find_gene_info, axis=1)

        new_results_path = os.path.abspath(os.path.join(gwas_results_file, '..'))
        fname = 'final_' +  association_name
        new_results_path = os.path.join(new_results_path, fname )
        gwas_results.to_csv(path_or_buf=new_results_path, sep='\t', index=False)
        description = "Genelist for GWAS results of trait " + association_name
         
        featureset_obj = self.build_featureset( new_results_path, genome_ref, description, workspace_name, association_name, prefix)
        
        return featureset_obj
