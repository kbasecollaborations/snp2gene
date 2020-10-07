/*
A KBase module: snp2gene
*/

module snp2gene {
    /*
        KBase style object reference X/Y/Z
            @id ws KBaseGwasData.Associations
    */
    typedef string association_ref;

    /*
        KBase style object reference X/Y/Z
            @id ws KBaseCollections.FeatureSet
    */
    typedef string featureset_ref;



    /*
        KBase style object reference X/Y/Z
            @id ws KBaseGenomes.Genome
    */
    typedef string genome_ref;




    /*
		A string representing a workspace name.
	*/
	typedef string workspace_name;

    /*
        A valid file path
    */
    typedef string file_path;

    typedef structure {
        file_path gwas_result_file;
        genome_ref genome_obj;
    } annotate_gwas_input;

    typedef structure {
        file_path snp_to_gene_list;
    } annotate_gwas_output;

    /*
        annotate_gwas_results:
            inputs:
                file path to gwas results
                genome object - with reference to GFF file

            outputs:
                TSV file represented by shock/handle ids and

    */
    funcdef annotate_gwas_results(annotate_gwas_input params) returns (annotate_gwas_output output) authentication required;


    
    typedef structure {
        string workspace_name;
        list <association_ref> associations;
        string p_value;
        string prefix;
    } annotate_gwas_app_input;

    typedef structure {
         string report_name;
         string report_ref;
        featureset_ref featureset_obj;
    } annotate_gwas_app_output;

    funcdef annotate_gwas_results_app(annotate_gwas_app_input params) returns (annotate_gwas_app_output output) authentication required;

};
