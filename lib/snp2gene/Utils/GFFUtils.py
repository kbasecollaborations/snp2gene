import tabix
import os
from pprint import pprint as pp

from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.GenomeFileUtilClient import GenomeFileUtil

class GFFUtils:
    def __init__(self, config):
        self.callback_url = config['SDK_CALLBACK_URL']
        self.shared_folder = config['scratch']

        self.dfu = DataFileUtil(self.callback_url)
        self.gfu = GenomeFileUtil(self.callback_url)

    def _get_dfu_obj(self, obj):
        return self.dfu.get_objects({'object_refs': [obj]})['data']

    def annotate_GWAS_results(self, genome_ref, association_ref):
        # gwas_results = self._get_dfu_obj(association_ref)
        # genome = self._get_dfu_obj(genome_ref)
        # gff_as_shock_id = self.gfu.export_genome_as_gff(genome_ref)

        wsid = self.dfu.ws_name_to_id('rmr:narrative_1552501344207')
        print(wsid)

        return 0
