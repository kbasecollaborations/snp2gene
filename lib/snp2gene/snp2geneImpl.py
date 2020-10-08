# -*- coding: utf-8 -*-
#BEGIN_HEADER
import logging
import os
import uuid
from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.WorkspaceClient import Workspace
from snp2gene.Utils.GFFUtils import GFFUtils
from snp2gene.Utils.GFFUtils2 import GFFUtils2
from installed_clients.KBaseReportClient import KBaseReport
#END_HEADER


class snp2gene:
    '''
    Module Name:
    snp2gene

    Module Description:
    A KBase module: snp2gene
    '''

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "0.0.1"
    GIT_URL = "git@github.com:kbasecollaborations/snp2gene.git"
    GIT_COMMIT_HASH = "8dd593e96c4b37fcf91a719181389e1b04c0bb4a"

    #BEGIN_CLASS_HEADER
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.config = config
        self.config['callback_url'] = os.environ['SDK_CALLBACK_URL']
        callback_url = self.config['callback_url']
        self.shared_folder = config['scratch']
        self.ws_url = config['workspace-url'] 
        self.wsc = Workspace(self.ws_url)
        self.kbr = KBaseReport(callback_url)
        logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
                            level=logging.INFO)
        #END_CONSTRUCTOR
        pass


    def annotate_gwas_results(self, ctx, params):
        """
        annotate_gwas_results:
        inputs:
            file path to gwas results
            genome object - with reference to GFF file
        outputs:
            TSV file represented by shock/handle ids and
        :param params: instance of type "annotate_gwas_input" -> structure:
           parameter "gwas_result_file" of type "file_path" (A valid file
           path), parameter "genome_obj" of type "genome_ref" (KBase style
           object reference X/Y/Z @id ws KBaseGenomes.Genome)
        :returns: instance of type "annotate_gwas_output" -> structure:
           parameter "snp_to_gene_list" of type "file_path" (A valid file
           path)
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN annotate_gwas_results

        gene_list = GFFUtils(self.config).annotate_GWAS_results(params['genome_obj'], params['gwas_result_file'])

        output = {'snp_to_gene_list': gene_list}

        #END annotate_gwas_results

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method annotate_gwas_results return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]

    def annotate_gwas_results_app(self, ctx, params):
        """
        :param params: instance of type "annotate_gwas_app_input" ->
           structure: parameter "associations" of list of type
           "association_ref" (KBase style object reference X/Y/Z @id ws
           KBaseGwasData.Associations), parameter "p_value" of String,
           parameter "prefix" of String
        :returns: instance of type "annotate_gwas_app_output" -> structure:
           parameter "report_name" of String, parameter "report_ref" of
           String, parameter "featureset_obj" of type "featureset_ref" (KBase
           style object reference X/Y/Z @id ws KBaseCollections.FeatureSet)
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN annotate_gwas_results_app
        # return the results
        print (params)
        #TODO: Hanlde cases where there are no significant SNPs
        #genome_ref = "47506/4/1"
        objects_created = []
        for association_ref in params['associations']:

            variation_ref = self.wsc.get_object_subset([{
                    'included': ['/variation_id'],
                    'ref': association_ref
                }])[0]['data']['variation_id']

            genome_ref = self.wsc.get_object_subset([{
                    'included': ['/genome_ref'],
                    'ref': variation_ref
                }])[0]['data']['genome_ref']

            featureset_obj  = GFFUtils2(self.config).annotate_GWAS_results(genome_ref, association_ref, params['workspace_name'], params['prefix'], params['p_value'])
            objects_created.append({'ref': featureset_obj,
                                        'description': 'FeatureSet'})
        # Build the new gff before doing anything
        
        # Download the workspace object for association one at a time
        # Filter SNPs for p-value, if no snps shows up, append this to warnings
        # Build the table structure needed for snp2gene
        # Run snp2gene algorithm and get final list.txt
        # Save as featureset. Find how to save featureset from genelist
        
        report_info = self.kbr.create_extended_report({
                'message': ' ',
                'objects_created': objects_created,
                'report_object_name': 'annotate_gwas_results_app_' + str(uuid.uuid4()),
                'workspace_name': params['workspace_name']
                })
        output = dict()
        output['report_name'] = report_info['name']
        output['report_ref'] = report_info['ref']
        print (output)

        #END annotate_gwas_results_app

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method annotate_gwas_results_app return value ' +
                             'output is not type dict as required.')
        return [output]
    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK",
                     'message': "",
                     'version': self.VERSION,
                     'git_url': self.GIT_URL,
                     'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
