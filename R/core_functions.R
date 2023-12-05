#' Retrive a desired PPI network from STRINGdb in a PPIpriority object
#' @description The build_ppi_network function is used to retrieve and
#' build a PPI network from the STRINGdb in a PPIpriority object. This
#' is the first step of the PPIpriority workflow, which allows to
#' perform random walk based target prioritization by calling
#' "prioritize_targets" on the output of this function.
#'
#' @import igraph STRINGdb
#'
#' @export
#' @usage build_ppi_network(version ="11.5", species = 9606,
#' score_threshold=400, input_directory=getwd())
#' @param version a character value, identifies the version of STRINGdb to use.
#' @param species a numeric value, that identifies the species of interest.
#' @param score_threshold a numeric value set to prune edges in the PPI network,
#' based on the level of confidence in the interaction. It should be a numeric
#' value between 0 and 1000. Values near to 1000 may cause problems with
#' prioritization, since will cause an intense pruning of the network.
#'
#' Note: score does not indicate the strength or the specificity of the
#' interaction. Instead, it is an indicator of confidence.
#'
#' @param input_directory Directory in which STRINGdb will download and
#' store the network information.
#' @return The output is an a PPIpriority_obj object (S4 class). The
#' object stores the STRINGdb connection and the PPI graph, that can
#' be accessed with the methods "get_STRINGdb_connection" and
#' "get_network", respectively.
#' The network is unweighted and the edges are already pruned for the
#' input value of score threshold. The vertices in this network have
#' also the attribute "preferred name" which is the HUGO nomeclature
#' for the STRING id in human network.
#'
#' Please refer to STRINGdb for version and species accessors.
#' If the function fails to download the required files because of having
#' reached download timeout, please extend timeout using "options(timeout=600)".
#' Removed the uncompleted downloaded files and call this function again.
#'
#' @examples
#' #Set timeout, not necessary when using the package
#' options(timeout=600)
#'
#' # Retrieve the PPI netowrk using STRINGdb
#' network_ppi_obj<- build_ppi_network(version ="11.5", species = 9606,
#' score_threshold=800, input_directory=getwd()) #Build a human PPI



build_ppi_network<-function(version ="11.5", species = 9606,
                            score_threshold=400, input_directory=getwd()){

  # input validation: predominantly evoked by STRINGdb
  check_input_threshold(score_threshold)

  # Retrieve the network from STRINGdb
  string_db_net<- STRINGdb::STRINGdb$new(version=version, species=species,
                                         score_threshold=score_threshold,
                                         input_directory=input_directory)

  # call the db (string_db_net alone does not work!)
  print(string_db_net)
  # get protein STRINGdb IDs
  protein_ids<- data.frame( STRING_id =names(V(string_db_net$graph)) )
  # add protein information or human, it retrieves mostly HUGO symbols
  ppi_net_info<-string_db_net$add_proteins_description(protein_ids)


  # Retrieve the COO graph format
  coo_network_format<- string_db_net$get_interactions(names(
    V(string_db_net$graph)))
  #removing duplicated edges
  coo_network_format<-unique(coo_network_format)
  colnames(coo_network_format)<- c("from","to","weight")

  # Build the graph using igraph
  ppi_network<-igraph::graph_from_data_frame(coo_network_format,
                                             directed = FALSE)

  # Add preferred gene name info label on the verteces (HUGO for human)
  index<-match(names(V(ppi_network)), ppi_net_info$STRING_id)
  ppi_network<-set.vertex.attribute(ppi_network, 'symbol',
                                    V(ppi_network),
                                    ppi_net_info$preferred_name[index])


  # Create PPIpriority obj (S4 class)
  PPIpriority_obj<-new("PPIpriority_obj",
                       STRINGdb_connection = string_db_net, graph = ppi_network)

  return(PPIpriority_obj)
}




#' Prioritize targets on PPI network via random walk from seed position(s).
#' @description Target prioritization based on Random Walk (RW) in the
#' selected STRINGdb PPI network based on seed potion(s).The function
#' requires to set as input seed genes and a STRINGdb PPI network,
#' which can be conveniently built using the function "build_ppi_network"
#' available in this package. Random Walk is based on the RWR function from
#' the RANKS package.
#'
#' @import igraph STRINGdb
#' @importFrom RANKS RWR
#' @importFrom dplyr filter select
#' @importFrom plyr join
#' @importFrom methods new
#' @export
#'
#' @usage prioritize_targets(PPIpriority_obj, seeds, candidates=NA, gamma= 0.6,
#'  tmax= 1000, eps=1e-10)
#' @param PPIpriority_obj A PPIpriority object.
#' @param seeds A vector of characters describing genes of interest
#' for seeding the random walk.
#' @param candidates A vector of characters describing genes which may
#' be associated to the phenotype of interest.
#' @param gamma A numeric value, restart parameter.
#' @param tmax A numeric value, maximum number of iterations.
#' @param eps A numeric value, maximum allowed difference between the
#' computed probabilities at the steady state.
#'
#'
#' @return data.frame of the target prioritization.
#' If the function fails to download the required files for downloading
#' matching annotations because of having reached download timeout,
#' please extend it using options(timeout=600)
#'
#' @references Wang et al. (2011) Network-based methods for human disease
#' gene prediction. Briefings in Functional Genomics 10(5):280–293.
#'
#' @examples
#' #Set extended timeout, not necessary while using the package
#' options(timeout=600)
#'
#' # Build a human PPI
#' network_ppi_obj<- build_ppi_network(version ="11.5", species = 9606,
#' score_threshold=800, input_directory=getwd())
#'
#' # Run the random walk using ERP44 and PRDX4 as seeds, print the output
#' prioritize_targets(network_ppi_obj, seeds=c("ERP44", "PRDX4"))


prioritize_targets<-function(PPIpriority_obj, seeds, candidates=NA,
                             gamma= 0.6, tmax= 1000, eps=1e-10){   #RW inputs


  #########################################
  # input valdiation
  check_PPIpriority_obj_input(PPIpriority_obj)  # PPIpriority_obj

  # test seeds/candidates inputs
  check_gene_names_input(seeds)
  if(sum(!is.na(candidates))>0){ check_gene_names_input(candidates)}

  # test numeric inputs for RWR
  check_numeric_input(gamma)
  check_numeric_input(tmax)
  check_numeric_input(eps)


  #########################################
  # retrieve STRINGdb ids for input gene seeds
  input_seed_genes<-data.frame(gene_names=seeds)
  conn<-get_STRINGdb_connection(PPIpriority_obj)    #Retrieve obj connection
  mapped_seeds_df<-conn$map(input_seed_genes,  "gene_names",
                            removeUnmappedRows = TRUE)


  check_mapping_output_seeds(mapped_seeds_df, seeds)

  # ... and for candidates, if any
  if(sum(!is.na(candidates))>0 ){
    input_candidates<-data.frame(gene_names=candidates)
    mapped_cand_df<-conn$map(input_candidates,  "gene_names",
                             removeUnmappedRows = TRUE)
    check_mapping_output_candidates(mapped_cand_df, candidates)
  }

  ########################################


  # extract adj matrix (just interactions, not weights)
  adj_matrix<-as_adjacency_matrix(get_network(PPIpriority_obj),
                                  type = "upper" , sparse = FALSE)

  # extract vertex names with STRING id and their associated preferred id
  vertex_names<-data.frame(STRING_id= V(get_network(PPIpriority_obj))$name,
                           preferred_id= V(get_network(PPIpriority_obj))$symbol)



  # check if seeds are in the network
  check_presence_seeds_in_graph(mapped_seeds_df, adj_matrix)

  if(sum(!is.na(candidates))>0 ){     # determine if cands are present
    # mapping df not 0 has been already checked
    check_presence_cand_in_graph(mapped_cand_df,adj_matrix)
  }



  # Random walk
  # use STRING_id mapped seeds
  res<-RANKS::RWR(adj_matrix, mapped_seeds_df$STRING_id,
                  gamma = gamma, tmax = tmax, eps = eps, norm = TRUE)

  #Remove seeds from output
  res_unsorted<-res$p[!names(res$p)%in%mapped_seeds_df$STRING_id]


  ranked_prob<-sort(res_unsorted[res_unsorted>0], decreasing = TRUE)
  rank_of_prob<-rank(-ranked_prob, ties.method="average")

  results<-data.frame(STRING_id=names(ranked_prob), rank=rank_of_prob,
                      prob_RW=ranked_prob)

  results<-plyr::join(results, vertex_names, by="STRING_id")
  results<-dplyr::select(results, c("STRING_id", "preferred_id",
                                    "rank", "prob_RW"))

  if(sum(!is.na(candidates))>0){
    check_presence_cand_in_res(mapped_cand_df, results)
    results<-filter(results, results$STRING_id %in% mapped_cand_df$STRING_id)
  }

  return(results)
}



#' PPI priority object
#'
#' @description PPI priority object. It stores STRINGdb and igraph objects.
#'
#' @slot STRINGdb_connection A STRINGdb object
#' @slot graph An igraph object
#' @importFrom methods setClass
#' @exportClass PPIpriority_obj
#' @examples
#' #Set timeout, not necessary while using the package
#' options(timeout=600)
#'
#' #Build a human PPI
#' network_ppi_obj<- build_ppi_network(version ="11.5", species = 9606,
#' score_threshold=800, input_directory=getwd())
#'
#' # Access network igraph
#' get_network(network_ppi_obj)
#'
#' # Access STRINGdb
#' get_STRINGdb_connection(network_ppi_obj)


setClass("PPIpriority_obj",
         slots=list(STRINGdb_connection="ANY", graph="ANY"))



#' Get PPI network from a PPIpriority_obj
#'
#' Accessor to the igraph object stored in a PPIpriority obj.
#' @importFrom methods setGeneric
#' @param object a PPIpriority object.
#' @return An igraph object
#' @examples
#' #Set extended timeout, not necessary while using the package
#' options(timeout=600)
#'
#' #Build human PPI
#' network_ppi_obj<- build_ppi_network(version ="11.5", species = 9606,
#' score_threshold=800, input_directory=getwd())
#'
#' # Access network igraph
#' get_network(network_ppi_obj)

setGeneric("get_network", function(object) {
  standardGeneric("get_network")
})




#' Get PPI network from a PPIpriority_obj
#'
#' Accessor to the igraph object stored in a PPIpriority obj.
#' @importFrom methods setMethod
#' @param object a PPIpriority object.
#' @docType methods
#' @rdname PPIpriority_obj-methods
#' @exportMethod get_network
#' @return An igraph object
#' @examples
#' #Set extended timeout, not necessary while using the package
#' options(timeout=600)
#'
#' #Build human PPI
#' network_ppi_obj<- build_ppi_network(version ="11.5", species = 9606,
#' score_threshold=800, input_directory=getwd())
#'
#' # Access network igraph
#' get_network(network_ppi_obj)


setMethod("get_network", signature("PPIpriority_obj"),
          definition = function(object){
            check_PPIpriority_obj_input(object)
            object@graph}
)




#' Get PPI network from a PPIpriority_obj
#'
#' Accessor to the STRINGdb obj stored in a PPIpriority obj.
#' @importFrom methods setGeneric
#' @param object a PPIpriority object.
#' @return A STRINGdb object.
#' @examples
#' #Set extended timeout, not necessary while using the package
#' options(timeout=600)
#'
#' #Build human PPI
#' network_ppi_obj<- build_ppi_network(version ="11.5", species = 9606,
#' score_threshold=800, input_directory=getwd())
#'
#' # Access STRINGdb
#' get_STRINGdb_connection(network_ppi_obj)

setGeneric("get_STRINGdb_connection", function(object) {
  standardGeneric("get_STRINGdb_connection")
})


#' Get PPI network from a PPIpriority_obj
#'
#' Accessor to the STRINGdb obj stored in a PPIpriority obj.
#' @importFrom methods setMethod
#' @docType methods
#' @rdname PPIpriority_obj-methods
#' @exportMethod get_STRINGdb_connection
#' @return A STRINGdb object
#' @examples
#' #Set extended timeout, not necessary while using the package
#' options(timeout=600)
#'
#' #Build a human PPI
#' network_ppi_obj<- build_ppi_network(version ="11.5", species =9606,
#' score_threshold=800, input_directory=getwd())
#'
#' # Access STRINGdb
#' get_STRINGdb_connection(network_ppi_obj)


setMethod("get_STRINGdb_connection", signature("PPIpriority_obj"),
          definition = function(object){
            check_PPIpriority_obj_input(object)
            object@STRINGdb_connection}
)
