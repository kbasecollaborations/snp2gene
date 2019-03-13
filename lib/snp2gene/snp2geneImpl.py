# -*- coding: utf-8 -*-
#BEGIN_HEADER
import logging
import os

from installed_clients.KBaseReportClient import KBaseReport

from snp2gene.Utils.GFFUtils import GFFUtils
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
    GIT_URL = ""
    GIT_COMMIT_HASH = ""

    #BEGIN_CLASS_HEADER
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.config = config
        self.config['callback_url'] = os.environ['SDK_CALLBACK_URL']
        self.shared_folder = config['scratch']
        logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
                            level=logging.INFO)
        #END_CONSTRUCTOR
        pass

    def annotate_gwas_results(self, ctx, params):
        """
        annotate_gwas_results:
        inputs:
            association object - with gwas results in a tuple
            genome object - with reference to GFF file
            workspace name
        outputs:
            TSV file represented by shock/handle ids and
        :param params: instance of type "annotate_gwas_input" -> structure:
           parameter "assoc_obj" of type "association_ref" (KBase style
           object reference X/Y/Z @id ws KBaseGwasData.Associations),
           parameter "genome_obj" of type "genome_ref" (KBase style object
           reference X/Y/Z @id ws KBaseGenomes.Genome), parameter
           "workspace_name" of type "workspace_name" (A string representing a
           workspace name.)
        :returns: instance of type "annotate_gwas_output" -> structure:
           parameter "snp_to_gene_list" of type "file_path" (A valid file
           path)
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN annotate_gwas_results

        gene_list = GFFUtils(self.config).annotate_GWAS_results(params['genome_obj'], params['assoc_obj'])

        output = {'snp_to_gene_list': '/path/to/snptogene/list'}
        #END annotate_gwas_results

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method annotate_gwas_results return value ' +
                             'output is not type dict as required.')
        # return the results
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
