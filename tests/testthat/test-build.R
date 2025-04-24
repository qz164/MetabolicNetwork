
#dir.create("inst/testdata", recursive = TRUE, showWarnings = FALSE)
library(testthat)

test_that("build_network returns proper structure and writes files", {
  data_path <- system.file("testdata", package = "MetabolicNetwork")

  #Merged_data <- read.delim(file.path(data_path, "mortpac_liver_metab_rna.txt"),
  #                          as.is = TRUE, check.names = FALSE, header = TRUE, na.strings = c("", " ", NA))
  #info <- read.delim(file.path(data_path, "info_integrate.txt"),
  #                   as.is = TRUE, check.names = FALSE, header = TRUE, na.strings = c("", " ", NA))

  #Merged_data = standardize_input(Merged_data, info)
  reactome.edge = read.delim(file.path(data_path, "Reactome_edges.txt"), as.is = T, check.names = F, header = T, na.strings = c("", " ", NA))
  reactome.node = read.delim(file.path(data_path, "Reactome_nodes.txt"), as.is = T, check.names = F, header = T, na.strings = c("", " ", NA))
  kegg.edge = read.delim(file.path(data_path, "KEGG_edges.txt"), as.is = T, check.names = F, header = T, na.strings = c("", " ", NA))
  kegg.node = read.delim(file.path(data_path, "KEGG_nodes.txt"), as.is = T, check.names = F, header = T, na.strings = c("", " ", NA))
  hmdb.edge = read.delim(file.path(data_path, "HMDB_edges.txt"), as.is = T, check.names = F, header = T, na.strings = c("", " ", NA))
  hmdb.node = read.delim(file.path(data_path, "HMDB_nodes.txt"), as.is = T, check.names = F, header = T, na.strings = c("", " ", NA))

  # Test with write_output = TRUE
  tmp_dir <- "_out_"
  dir.create(tmp_dir)

  #result = build_kegg(
  #  Merged_data = Merged_data,
  #  info = info,
  #  species = "rno",       # or "mmu"/"hsa"
  #  gene_id = "Gene name",   # column name for genes
  #  met_id = "KEGG",      # column name for metabolites
  #  size = "MetMN_Expand",
  #  write_output = T,
  #  output_dir = tmp_dir
  #)

  #result = build_hmdb(
  #  Merged_data = Merged_data,
  #  info = info,
  #  species = "rno",       # or "mmu"/"hsa"
  #  gene_id = "Gene name",   # column name for genes
  #  met_id = "HMDB lip",      # column name for metabolites
  #  size = "Diff_Expand",
  #  write_output = T,
  #  output_dir = tmp_dir
  #)
  result = merge_network(
    reactome.edge,
    reactome.node,
    kegg.edge,
    kegg.node,
    hmdb.edge,
    hmdb.node,
    write_output = T,
    output_dir = tmp_dir
  )


  files <- list.files(tmp_dir)
  expect_true(all(c("Merged_nodes.txt",
                    "Merged_edges.txt",
                    "Merged_edges_solid.txt") %in% files))
})
