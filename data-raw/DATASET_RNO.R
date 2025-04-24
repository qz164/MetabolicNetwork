## code to prepare data resources

reactome.reaction.rno = read.delim("data-raw/Reactome_reaction_RNO.txt", as.is = T, check.names = F, header = T, na.strings = c("", " ", NA))
ont_chebi = read.csv("data-raw/chebi_ontology.csv", as.is = T, check.names = F, header = T)
chebi = read.delim("data-raw/ChEBI_data.tsv", as.is = T, check.names = F, header = T, quote = "")
map = read.csv("data-raw/compoundID_map.csv", as.is = T, check.names = F, header = T, na.strings = c("", " ", NA))

kegg.reaction.rno = read.delim("data-raw/KEGG_reaction_RNO.txt", as.is = T, check.names = F, header = T, na.strings = c("", " ", NA))
HMDB.reaction.rno = read.delim("data-raw/HMDB_protein_association_Rat.txt", as.is = T, check.names = F, header = T, na.strings = c(" ", "", NA))

lipmap = readxl::read_excel("data-raw/HMDB-map-lipids.xlsx")
usethis::use_data(reactome.reaction.rno,
                  kegg.reaction.rno,
                  HMDB.reaction.rno,
                  ont_chebi, chebi, map,
                  lipmap,
                  internal = TRUE,
                  overwrite = TRUE)



