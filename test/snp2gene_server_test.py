# -*- coding: utf-8 -*-
import os
import time
import unittest
import shutil
from configparser import ConfigParser

from snp2gene.snp2geneImpl import snp2gene
from snp2gene.snp2geneServer import MethodContext
from snp2gene.authclient import KBaseAuth as _KBaseAuth

from installed_clients.WorkspaceClient import Workspace


class snp2geneTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        token = os.environ.get('KB_AUTH_TOKEN', None)
        config_file = os.environ.get('KB_DEPLOYMENT_CONFIG', None)
        cls.cfg = {}
        config = ConfigParser()
        config.read(config_file)
        for nameval in config.items('snp2gene'):
            cls.cfg[nameval[0]] = nameval[1]
        # Getting username from Auth profile for token
        authServiceUrl = cls.cfg['auth-service-url']
        auth_client = _KBaseAuth(authServiceUrl)
        user_id = auth_client.get_user(token)
        # WARNING: don't call any logging methods on the context object,
        # it'll result in a NoneType error
        cls.ctx = MethodContext(None)
        cls.ctx.update({'token': token,
                        'user_id': user_id,
                        'provenance': [
                            {'service': 'snp2gene',
                             'method': 'please_never_use_it_in_production',
                             'method_params': []
                             }],
                        'authenticated': 1})
        cls.wsURL = cls.cfg['workspace-url']
        cls.wsClient = Workspace(cls.wsURL)
        cls.serviceImpl = snp2gene(cls.cfg)
        cls.scratch = cls.cfg['scratch']
        cls.callback_url = os.environ['SDK_CALLBACK_URL']
        suffix = int(time.time() * 1000)
        cls.wsName = "test_ContigFilter_" + str(suffix)
        ret = cls.wsClient.create_workspace({'workspace': cls.wsName})  # noqa

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')

    # NOTE: According to Python unittest naming rules test method names should start from 'test'. # noqa
    def test_your_method(self):
        # Prepare test objects in workspace if needed using
        # self.getWsClient().save_objects({'workspace': self.getWsName(),
        #                                  'objects': []})
        #
        # Run your method by
        # ret = self.getImpl().your_method(self.getContext(), parameters...)
        #
        # Check returned data with
        # self.assertEqual(ret[...], ...) or other unittest methods

        result_file = shutil.copyfile('/kb/module/test/sample_data/snpdata-poplar.tsv', os.path.join(self.scratch, 'snpdata.tsv'))

        ret = self.serviceImpl.annotate_gwas_results(self.ctx, {
            'workspace_name': 'rmr:narrative_1552501344207',
            # Populus v3.1
            'genome_obj': '26587/12/1',
            'gwas_result_file': os.path.join(self.scratch, 'snpdata.tsv')

            # A. Thaliana TAIR10
            #'genome_obj': '25404/2/1',
            #'gwas_result_file': '/kb/module/test/sample_data/sample_gwas_results.tsv'
        })
