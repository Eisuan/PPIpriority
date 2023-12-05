
#' Check PPI priority object input
#' @description Input validation function: called during the input validation
#' step of "prioritize_targets".
#' @import methods
#' @import STRINGdb
#' @import igraph
#' @param obj_input A PPIpriority_obj
#' @return An appropriate error message


check_PPIpriority_obj_input<-function(obj_input){
  if(!isS4(obj_input)){
    stop("The input is not an S4 class object. Please, use the
    \"build_ppi_network\" from this package to build the correct input for the
    function \"target_prioritization\" ")

  }else if(!is(obj_input,"PPIpriority_obj") ){
    stop("The S4 class input is not a PPIpriority_obj. Please, use
    the \"build_ppi_network\" from this package to build the correct
    input for the function \"target_prioritization\"")

  }else if(!igraph::is.igraph(obj_input@graph)){
    stop("The input object does not store the expected PPI igraph network.\n
    Please, use the \"build_ppi_network\" from this package to build the
    correct input for the function \"target_prioritization\" ")
  }

  else if(!is(obj_input@STRINGdb_connection, "STRINGdb")){
    stop("The input object does not store the expected STRINGdb object.\n
    Please, use the \"build_ppi_network\" from this package to build the
    correct input for the function \"target_prioritization\" ")
  }

}



#' Check gene names input
#' @description Input validation function: called during the input validation
#' step of "prioritize_targets".
#'
#' @param input_gene_names A vector of gene names
#' @return An appropriate error message

check_gene_names_input<-function(input_gene_names){
  if(!is.character(input_gene_names)){
    stop("The input gene list should be a character vector")
  }
}


#' Warning for high score_threshold
#' @description Input validation function: called during the input validation
#' step of "prioritize_targets".
#'
#' @param threshold A vector of gene names
#' @return An appropriate error message

check_input_threshold<-function(threshold){
  if(threshold>900){
    warning("Input score_threshold is higher than 900. High score threshold
    may retrieve no networks or cause downstrem problems in the analysis because
    of excessive pruning of the PPI netowork. For Homo Sapiens PPI, values equal
    to or higher than 900 start to prune eccesively the network. Values over
    1000 are not allowed.")
  }
}


#' Check numeric inputs
#' @description Input validation function: called during the input validation
#' step of "prioritize_targets"to check RWR inputs.
#'
#' @param x an input parameter to check
#' @return An appropriate error message

check_numeric_input<-function(x){
  if(!is.numeric(x)){
    stop("The input is not numeric. Please, provide a suitable input for RW parameters")}
}




#' Check mapping output for seeds
#' @description Input validation function: called after mapping of input gene
#' names (seed) in "prioritize_targets". Prevents
#' "prioritize_targets" to continue if input genes can not be mapped to the
#' STRING ids of the PPI network. It also informs the user which are the genes
#' that are excluded from the analysis because of unsuccessful mapping.
#' The warning is raised by STRINGdb map.
#' @param df A data.frame containing the columns: "gene_names" , "STRING_id"
#' @param input_gene_names A vector of gene names
#' @return An appropriate error message or a complementary message to help
#' the user understanding the unmapped gene names.


check_mapping_output_seeds<-function(df, input_gene_names){
  if(nrow(df)==0){
    stop("STRINGdb failed to map all your seed gene inputs. This package utilizes
    STRINGdb map to enstablish the relationship between STRING id and the
    input gene list vector. Please, check the integrity of your list or refer
    to STRINGdb map method documentation for further information.")

  }else if(!nrow(df)==length(input_gene_names)){
    message("\nThe following genes have not been mapped:\n",
    paste(input_gene_names[!input_gene_names%in%df$gene_names], collapse=' '),
    "\nThus, they will not be considered during the analysis as
    seeds respectively.")
  }
}


#' Check mapping output for candidates
#' @description Input validation function: called after mapping of input gene
#' names (candidate) in "prioritize_targets". Prevents
#' "prioritize_targets" to continue if input genes can not be mapped to the
#' STRING ids of the PPI network. It also informs the user which are the genes
#' that are excluded from the analysis because of unsuccessful mapping.
#' The warning is raised by STRINGdb map.
#' @param df A data.frame containing the columns: "gene_names" , "STRING_id"
#' @param input_gene_names A vector of gene names
#' @return An appropriate error message or a complementary message to help
#' the user understanding which are the unmapped gene names.


check_mapping_output_candidates<-function(df, input_gene_names){
  if(nrow(df)==0){
    stop("STRINGdb failed to map all your candidate gene inputs. This package
    utilizes STRINGdb map to enstablish the relationship between STRING id and
    the input gene list vector. Please, check the integrity of your list or
    refer to STRINGdb map method documentation for further information.")

  }else if(!nrow(df)==length(input_gene_names)){
    message("\nThe following genes have not been mapped:\n",
            paste(input_gene_names[!input_gene_names%in%df$gene_names],
                  collapse=' '),
            "\nThus, they will not be considered during the analysis as
    candidates respectively.")
  }
}






#' Check presence for seeds in the graph
#' @description Internal validation function:  Prevents
#' "prioritize_targets" to continue if mapped input gene id can not be found
#' in the PPI network.
#' @param seed_df A data.frame of mapped seed ids from STRINGdb
#' containing the columns: "gene_names" , "STRING_id"
#' @param adj_matrix A matrix, adj. matrix from igraph obj in ppi_priority obj
#' @return An appropriate error message or a complementary message to help
#' the user understanding which are the missing seeds.


check_presence_seeds_in_graph<-function(seed_df, adj_matrix){
  if(sum(seed_df$STRING_id%in%colnames(adj_matrix)) != nrow(seed_df) ){

    seeds_missing<-seed_df$gene_names[!seed_df$STRING_id%in%colnames(adj_matrix)]
    message("\nThe following seeds have not been found in the graph:\n",
            paste(seeds_missing,collapse=' ') )

    stop("Failed to find your seeds in the graph. This is may be due to very
    high score thresholding values (near to 1000). Please, build a new PPI
    object with a lower score_threshold parameter (e.g. 400) OR remove seeds
    not found in the graph.")
  }
}



#' Check presence for candidates in the graph
#' @description Internal validation function:  Prevents
#' "prioritize_targets" to continue if mapped input gene id can not be found
#' in the PPI network.
#' @param can_df A data.frame of mapped candidate ids from STRINGdb
#' containing the columns: "gene_names" , "STRING_id"
#' @param adj_matrix A matrix, adj. matrix from igraph obj in ppi_priority obj
#' @return An appropriate error message or a complementary message


check_presence_cand_in_graph<-function(can_df, adj_matrix){
  if(sum(can_df$STRING_id%in%colnames(adj_matrix)) != nrow(can_df) ){

    can_missing<-can_df$gene_names[!can_df$STRING_id%in%colnames(adj_matrix)]
    message("\nThe following candidates have not been found in the graph:\n",
            paste(can_missing,collapse=' ') )

    stop("Failed to find your candidates in the graph. This is may be due to
    selection of very high score thresholding values (near to 1000) while
    creating the ppi_priority object. Please, build a new PPI object with
    a lower score_threshold parameter (e.g. 400) OR remove candidates
    genes not found in the graph.")

  }
}



#' Check presence for candidates in the RW results
#' @description Internal validation function:  warns that candidates are not
#' in RW results.
#' @param can_df A data.frame of mapped candidate ids from STRINGdb
#' containing the columns: "gene_names" , "STRING_id"
#' @param results A data.frame, internal result of prioritize_targets.
#' @return An appropriate warning


check_presence_cand_in_res<-function(can_df, results){
  if(sum(can_df$STRING_id %in% results$STRING_id ) != nrow(can_df) ){

    can_missing<-can_df$gene_names[!can_df$STRING_id%in%results$STRING_id]
    warning("\nThe following candidates have not scored significant RW results:\n",
    paste(can_missing,collapse=' '),
    "\nIf no candidate genes are showed but you are still interested
    in the other results, please perform the analysis without candidates")

  }
}


