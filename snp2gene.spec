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
        association_ref assoc_obj;
        genome_ref genome_obj;
        workspace_name workspace_name;
    } annotate_gwas_input;

    typedef structure {
        file_path snp_to_gene_list;
    } annotate_gwas_output;

    /*
        annotate_gwas_results:
            inputs:
                association object - with gwas results in a tuple
                genome object - with reference to GFF file
                workspace name

            outputs:
                TSV file represented by shock/handle ids and

    */
    funcdef annotate_gwas_results(annotate_gwas_input params) returns (annotate_gwas_output output) authentication required;

};
