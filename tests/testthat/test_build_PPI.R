test_that("Evoke input errors/messages", {
  #extend timeout
  options(timeout=600)

  ## Total number of tests 11

  # Load a dataset with default parameters
  ppi_net<-build_ppi_network()

  #### Testing build_ppi_network and STRINGdb error invocation ####

  ### 1. version
  # 1.1 Numeric input instead of string (STRINGdb error)
  expect_error(build_ppi_network(version =11.5, species = 9606,
                                 score_threshold=800, input_directory=getwd()))

  # 1.2 Call a non existing version (STRINGdb error)
  expect_error(build_ppi_network(version ="23d2", species = 9606,
                                 score_threshold=800, input_directory=getwd()))

  ## 2. species
  # 2.1 String input instead of numeric (STRINGdb error)
  expect_error(build_ppi_network(version ="11.5", species = "9606",
                                 score_threshold=800, input_directory=getwd()))

  # 2.2 Call a non existing version  (STRINGdb error)
  expect_error(build_ppi_network(version ="11.5", species = 1164864846,
                                 score_threshold=800, input_directory=getwd()))


  ## 3. score_threshold
  # 3.1 Call a score version as a string (STRINGdb error)
  expect_error(build_ppi_network(version ="11.5", species = 9606,
                                score_threshold="800",input_directory=getwd()))

  expect_warning(build_ppi_network(version ="11.5", species = 9606,
                                 score_threshold=950,input_directory=getwd()))

  ## 4 input_directory
  # 4.1 Call a directory as a number
  expect_error(build_ppi_network(version ="11.5", species = 9606,
                                 score_threshold=800, input_directory=5874))



  #### Testing prioritize_targets ####

  # Arguments: PPIpriority_obj, seeds, candidates, gamma, tmax,
  # eps, print_plots

  setClass("fake_class", slots=list(name="character", surname="character"))
  fakes4<- new("fake_class", name = "Petro", surname = "Frittata")

  faulty_creation_ppi_priority<-new("PPIpriority_obj",
                           STRINGdb_connection = "Petro", graph= "Frittata")

  # Creating broken versions of the ppi_priority_obj
  ppi_priority_STRINGdb_breaking<-ppi_net
  ppi_priority_STRINGdb_breaking@STRINGdb_connection<-"bug"

  ppi_priority_graph_breaking<-ppi_net
  ppi_priority_graph_breaking@graph<-"bug"

  seed=c("ERP44", "LMAN1")
  seed_faulty_mixed_type=c("ERP44", 35423)
  seed_faulty_wrong_type=c(6985, 35423)

  # 0. No inputs
  expect_error(res_df<-prioritize_targets())


  ## 1. PPIpriority obj
  # 1.1 Testing string and numeric outputs
  expect_error(prioritize_targets("string", seed ))
  expect_error(prioritize_targets(22, seed ))

  # 1.2 Testing S4 class: S4 class not ppi_priority obj
  expect_error(prioritize_targets(fakes4, seed ))

  # 1.3 Testing S4 class: S4 class ppi_priority with faulty creation
  expect_error(prioritize_targets(faulty_ppi_priority, seed))

  # 1.4 Testing S4 class: S4 class ppi_priority with broken STRINGdb slot
  expect_error(prioritize_targets(ppi_priority_STRINGdb_breaking,seed))

  # 1.5 Testing S4 class: S4 class ppi_priority with broken graph slot
  expect_error(prioritize_targets(ppi_priority_graph_breaking,seed))



  ## 2. Seeds
  ppi_net_pruned<- build_ppi_network(score_threshold=999) # for 3.3 test

  # 2.1 test 1 non mappable gene over 2 for seeds, expect message
  expect_message(prioritize_targets(ppi_net,seeds = c("ERP44", "NO_GENE")) )

  # 2.2 test non existing genes
  expect_error(prioritize_targets(ppi_net,seeds = c("NO_GENE1", "NO_GENE2")) )

  # 3.3 test absence of mapped genes in the final graph.
  # to to this, create a strongly pruned graph (which is not advisable)
  expect_error(prioritize_targets(ppi_net_pruned,
                                  seeds = c("ERP44", "MXD1", "DEDD",
                                            "LY9", "USF1")) )

  ## 3. Candidates

  # 2.1 test 1 non mappable gene over 2 for candidates, expect message
  expect_message(prioritize_targets(ppi_net,seeds = c("ERP44"),
                                    candidates = c("ERAP1", "NOGENE2")))

  # 2.2 test non existing genes
  expect_error(prioritize_targets(ppi_net,seeds = c("ERP44"),
                                  candidates = c("NOGENE1", "NOGENE2") ) )

  # 3.3 test absence of mapped genes in the final graph.
  # to to this, create a strongly pruned graph (which is not advisable)
  expect_error(prioritize_targets(ppi_net_pruned,
                                  seeds = c("TP53"),  # p53 is in the graph
                                  candidates = c("ERP44", "MXD1", "DEDD",
                                                 "LY9", "USF1")) )

  # check gamma, tmax and eps inputs
  expect_error(prioritize_targets(ppi_net_pruned, seeds = c("TP53"),
                                  gamma="0.6"))
  expect_error(prioritize_targets(ppi_net_pruned, seeds = c("TP53"),
                                  tmax="0.6"))
  expect_error(prioritize_targets(ppi_net_pruned, seeds = c("TP53"),
                                  eps="0.6"))


  # 4 test ppi priority methods
  # 4.1 Testing S4 class: S4 class ppi_priority with faulty creation
  expect_error(get_network(faulty_ppi_priority, seed))

  # 4.2 Testing S4 class: S4 class ppi_priority with broken STRINGdb slot
  expect_error(get_network(ppi_priority_STRINGdb_breaking))

  # 4.3 Testing S4 class: S4 class ppi_priority with broken graph slot
  expect_error(get_network(ppi_priority_graph_breaking))

  # 4.4 Testing S4 class: S4 class ppi_priority with faulty creation
  expect_error(get_STRINGdb_connection(faulty_ppi_priority, seed))

  # 4.5 Testing S4 class: S4 class ppi_priority with broken STRINGdb slot
  expect_error(get_STRINGdb_connection(ppi_priority_STRINGdb_breaking))

  # 4.6 Testing S4 class: S4 class ppi_priority with broken graph slot
  expect_error(get_STRINGdb_connection(ppi_priority_graph_breaking))




})
