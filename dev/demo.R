devtools::document()
devtools::install(force = TRUE)

library(MetabolicNetwork)
source("data-raw/DATASET_RNO.R")

data_path = system.file("testdata", package = "MetabolicNetwork")
Merged_data = read.delim(file.path(data_path, "mortpac_liver_metab_rna.txt"),
                         as.is = T, check.names = F, header = T, na.strings = c("", " ", NA))
info = read.delim(file.path(data_path, "info_integrate.txt"),
                  as.is = T, check.names = F, header = T, na.strings = c("", " ", NA))

# Standardize column names
Merged_data = standardize_input(Merged_data, info)

# Build the network
result = build_reactome(
  Merged_data = Merged_data,
  info = info,
  species = "rno",       # or "mmu"/"hsa"
  gene_id = "Uniprot",   # column name for genes
  met_id = "ChEBI",      # column name for metabolites
  size = "MetMN_Expand",
  write_output = T,
  output_dir = "outputs"
)

result = build_kegg(
  Merged_data = Merged_data,
  info = info,
  species = "rno",       # or "mmu"/"hsa"
  gene_id = "Gene name",   # column name for genes
  met_id = "KEGG",      # column name for metabolites
  size = "MetMN_Expand",
  write_output = T,
  output_dir = "outputs"
)

result = build_hmdb(
  Merged_data = Merged_data,
  info = info,
  species = "rno",       # or "mmu"/"hsa"
  gene_id = "Gene name",   # column name for genes
  met_id = "HMDB lip",      # column name for metabolites
  size = "Diff_Expand",
  write_output = T,
  output_dir = "outputs"
)


reactome.edge = read.delim("outputs/Reactome_edges.txt", as.is = T, check.names = F, header = T, na.strings = c("", " ", NA))
reactome.node = read.delim("outputs/Reactome_nodes.txt", as.is = T, check.names = F, header = T, na.strings = c("", " ", NA))
kegg.edge = read.delim("outputs/KEGG_edges.txt", as.is = T, check.names = F, header = T, na.strings = c("", " ", NA))
kegg.node = read.delim("outputs/KEGG_nodes.txt", as.is = T, check.names = F, header = T, na.strings = c("", " ", NA))
hmdb.edge = read.delim("outputs/HMDB_edges.txt", as.is = T, check.names = F, header = T, na.strings = c("", " ", NA))
hmdb.node = read.delim("outputs/HMDB_nodes.txt", as.is = T, check.names = F, header = T, na.strings = c("", " ", NA))

result = merge_network(
  reactome.edge,
  reactome.node,
  kegg.edge,
  kegg.node,
  hmdb.edge,
  hmdb.node,
  write_output = T,
  output_dir = "outputs"
)






