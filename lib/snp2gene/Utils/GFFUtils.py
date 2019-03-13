import tabix
import os
from pprint import pprint as pp

from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.GenomeFileUtilClient import GenomeFileUtil

class GFFUtils:
    def __init__(self, config):
        self.callback_url = config['callback_url']
        self.shared_folder = config['scratch']
        self.GFF_dir = os.path.join(self.shared_folder, 'GFF')

        if not os.path.isdir(self.GFF_dir):
            os.mkdir(self.GFF_dir)

        self.dfu = DataFileUtil(self.callback_url)
        self.gfu = GenomeFileUtil(self.callback_url)

    def annotate_GWAS_results(self, genome_ref, association_ref):
        gff = self.gfu.genome_to_gff({
            'genome_ref': genome_ref,
            'target_dir': self.GFF_dir
        })['file_path']

        tb = tabix.open(gff)

        assoc_results = self.dfu.get_objects({'object_refs': [association_ref]})['data'][0]['data']['association_details']

        pp(tb)

        return 0
