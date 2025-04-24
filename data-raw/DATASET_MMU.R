## code to prepare data resources

reactome.reaction.mmu = read.delim("data-raw/Reactome_reaction_MMU.txt", as.is = T, check.names = F, header = T, na.strings = c("", " ", NA))
ont_chebi = read.csv("data-raw/chebi_ontology.csv", as.is = T, check.names = F, header = T)
chebi = read.delim("data-raw/ChEBI_data.tsv", as.is = T, check.names = F, header = T, quote = "")
map = read.csv("data-raw/compoundID_map.csv", as.is = T, check.names = F, header = T, na.strings = c("", " ", NA))

kegg.reaction.mmu = read.delim("data-raw/KEGG_reaction_MMU.txt", as.is = T, check.names = F, header = T, na.strings = c("", " ", NA))
HMDB.reaction.mmu = read.delim("data-raw/HMDB_protein_association_Mouse.txt", as.is = T, check.names = F, header = T, na.strings = c(" ", "", NA))

lipmap = readxl::read_excel("data-raw/HMDB-map-lipids.xlsx")

usethis::use_data(reactome.reaction.mmu,
                  kegg.reaction.mmu,
                  HMDB.reaction.mmu,
                  ont_chebi, chebi, map,
                  lipmap,
                  internal = TRUE,
                  overwrite = TRUE)



