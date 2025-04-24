rm(list = ls()[])
library(dplyr)
library(tidyr)

ChEBI = read.delim("ChEBI_data.tsv", as.is = T, check.names = F, header = T, na.strings = c("", " ", NA), quote = "")
ChEBI$ChEBI_ID = gsub("CHEBI:", "", ChEBI$ChEBI_ID)
HMDB = read.csv("../HMDB/hmdb_2024.csv", as.is = T, check.names = F, header = T, na.strings = c("", " ", NA))
Lmaps = read.csv("LIPIDMAPS.csv", as.is = T, check.names = F, header = T, na.strings = c("", " ", NA))
kegg = read.delim("Kegg_compound.txt", as.is = T, check.names = F, header = F, na.strings = c("", " ", NA))
colnames(kegg) = c("ID", "NameAll")
kegg$Name = unlist(sapply(kegg$NameAll, function(x)strsplit(x, ";")[[1]][1]))


# KEGG - CHEBI id conversion map
chebi_kegg = read.delim("../kegg/chebi_kegg.txt", as.is = T, check.names = F, header = F)
colnames(chebi_kegg) = c("ChEBI", "KEGG")
chebi_kegg$ChEBI = gsub("chebi:", "", chebi_kegg$ChEBI)
chebi_kegg$KEGG = gsub("cpd:", "", chebi_kegg$KEGG)
kegg = kegg %>% full_join(chebi_kegg, by = c("ID" = "KEGG")) %>% distinct()


# Get a list of unique InChikey and add identifiers
table(ChEBI$InChIKey == "None") # 4710 /31702
table(is.na(HMDB$InChIKey)) # 21 / 217920
table(is.na(Lmaps$INCHI_KEY)) # no NA
test = ChEBI[ChEBI$InChIKey == "None",] # some of those are glycans

InChikey = unique(c(ChEBI$InChIKey, HMDB$InChIKey, Lmaps$INCHI_KEY))
InChikey = InChikey[!is.na(InChikey) & InChikey != "None"]
map = as.data.frame(InChikey)

map = map %>% left_join(HMDB[,c("accession", "name", "kegg", "chebi_id", "pubchem", "InChIKey")], 
                              by = c("InChikey" = "InChIKey"))
map = map %>% left_join(Lmaps[,c("LM_ID", "NAME", "INCHI_KEY", "PUBCHEM_CID", "CHEBI_ID")],
                        by = c("InChikey" = "INCHI_KEY"))
map = map %>% left_join(ChEBI[,c("ChEBI_ID", "ChEBI_Name", "HMDB Database Links", "KEGG_Links", "InChIKey")],
                        by = c("InChikey" = "InChIKey"))
map$CHEBI_ID = as.character(map$CHEBI_ID)
map$chebi_id = as.character(map$chebi_id)

# Add all KEGG id 
kegg = kegg[,-2]
colnames(kegg) = c("KEGG", "KEGG_Name", "ChEBI")
map = map %>% full_join(kegg, by = c("kegg" = "KEGG")) 

# Merge the identifiers with multiple columns
map = map %>%
  mutate(row_id = row_number()) %>%
  mutate(across(everything(), as.character)) %>%
  pivot_longer(cols = c(accession, `HMDB Database Links`), 
               names_to = "source", values_to = "HMDB") %>%
  group_by(row_id) %>%
  filter(!is.na(HMDB) | all(is.na(HMDB))) %>%
  ungroup() %>%
  select(-source, -row_id)
map = map[!duplicated(map),]

map = map %>% 
  mutate(row_id = row_number()) %>%
  mutate(across(everything(), as.character)) %>%
  pivot_longer(cols = c(CHEBI_ID, ChEBI_ID, chebi_id, ChEBI), 
               names_to = "source", values_to = "ChEBI") %>%
  group_by(row_id) %>%
  filter(!is.na(ChEBI) | all(is.na(ChEBI))) %>%
  ungroup() %>%
  select(-source, -row_id)
map = map[!duplicated(map),]

map = map %>% 
  mutate(row_id = row_number()) %>%
  mutate(across(everything(), as.character)) %>%
  pivot_longer(cols = c(kegg, KEGG_Links), 
               names_to = "source", values_to = "KEGG") %>%
  group_by(row_id) %>%
  filter(!is.na(KEGG) | all(is.na(KEGG))) %>%
  ungroup() %>%
  select(-source, -row_id)
map = map[!duplicated(map),]

map = map %>% 
  mutate(row_id = row_number()) %>%
  mutate(across(everything(), as.character)) %>%
  pivot_longer(cols = c(pubchem, PUBCHEM_CID), 
               names_to = "source", values_to = "PubChem") %>%
  group_by(row_id) %>%
  filter(!is.na(PubChem) | all(is.na(PubChem))) %>%
  ungroup() %>%
  select(-source, -row_id)
map = map[!duplicated(map),]


########################### Clean up names
colnames(map)[colnames(map) == "name"] = "HMDB_Name"
colnames(map)[colnames(map) == "NAME"] = "LM_Name"

map = map[,c("InChikey", "HMDB_Name", "LM_Name", "ChEBI_Name", "KEGG_Name",
             "HMDB", "ChEBI", "KEGG", "LM_ID", "PubChem")]
# An unique InChiKey usually can map to multiple ChEBI IDs
# An unique HMDB id sometimes map to multiple InChiKey (with first 14 characters the same)

map$Name = NA

# Group by HMDB id --> fill the HMDB id name 
map = map %>% group_by(HMDB) %>%
  mutate(Name = HMDB_Name) %>% 
  fill(Name, .direction = "downup")

# for those without HMDB id --> go for the shortest non-NA names 
for (i in 1:nrow(map)){
  if (is.na(map$Name[i])){
    tmp = as.character(map[i,grep("Name", colnames(map))])
    tmp = tmp[!is.na(tmp)]
    tmp = tmp[which.min(nchar(tmp))]
    map$Name[i] = ifelse(length(tmp) > 0, tmp, NA)
  }
}

map = map[,c("InChikey", "Name", "HMDB_Name", "LM_Name", "ChEBI_Name", "KEGG_Name", 
             "HMDB", "ChEBI", "KEGG", "LM_ID", "PubChem")]

########################### Fill up NAs
map = map %>%
  group_by(HMDB) %>%
  mutate(
    KEGG = ifelse(is.na(HMDB), KEGG, zoo::na.locf(KEGG, na.rm = FALSE)),
    ChEBI = ifelse(is.na(HMDB), ChEBI, zoo::na.locf(ChEBI, na.rm = FALSE))
  ) %>%
  ungroup() %>%
  group_by(KEGG) %>%
  mutate(
    HMDB = ifelse(is.na(KEGG), HMDB, zoo::na.locf(HMDB, na.rm = FALSE)),
    ChEBI = ifelse(is.na(KEGG), ChEBI, zoo::na.locf(ChEBI, na.rm = FALSE))
  ) %>%
  ungroup() %>%
  group_by(ChEBI) %>%
  mutate(
    KEGG = ifelse(is.na(ChEBI), KEGG, zoo::na.locf(KEGG, na.rm = FALSE)),
    HMDB = ifelse(is.na(ChEBI), HMDB, zoo::na.locf(HMDB, na.rm = FALSE)),
  ) %>%
  ungroup()

m = match(map$HMDB, HMDB$accession)
map$HMDB_Name = HMDB$name[m]

m = match(map$ChEBI, ChEBI$ChEBI_ID)
map$ChEBI_Name = ChEBI$ChEBI_Name[m]

m = match(map$KEGG, kegg$KEGG)
map$KEGG_Name = kegg$KEGG_Name[m]


###########################  Final Check --> Group by Name (all lower case) --> fill all IDs
map = map %>%
  mutate(Name_lower = tolower(Name)) %>%  
  group_by(Name_lower) %>%
  mutate(
    KEGG = ifelse(is.na(Name), KEGG, zoo::na.locf(KEGG, na.rm = FALSE)),
    ChEBI = ifelse(is.na(Name), ChEBI, zoo::na.locf(ChEBI, na.rm = FALSE)),
    HMDB = ifelse(is.na(Name), HMDB, zoo::na.locf(HMDB, na.rm = FALSE)),
    InChikey = ifelse(is.na(Name), InChikey, zoo::na.locf(InChikey, na.rm = FALSE)),
  ) %>%
  ungroup() %>%
  select(-Name_lower) %>%  
  distinct(KEGG, HMDB, ChEBI, InChikey, LM_ID, .keep_all = TRUE)




################################################################## 
# Another map to convert compounds + change back to the original compound

chebi.charge = ChEBI[ChEBI$Charge != 0,]
chebi.charge$Charge = ifelse(chebi.charge$Charge < 0, paste0(abs(chebi.charge$Charge), "-"), paste0(chebi.charge$Charge, "+"))
chebi.charge$Charge = paste0("(", chebi.charge$Charge, ")")

chebi.charge$ChEBI_Name_Convert = mapply(function(a, b) gsub(a, "", b, fixed = TRUE), chebi.charge$Charge, chebi.charge$ChEBI_Name)
chebi.charge = chebi.charge[chebi.charge$ChEBI_Name != chebi.charge$ChEBI_Name_Convert,]

ii = intersect(chebi.charge$ChEBI_Name_Convert, map$ChEBI_Name)
ii = ii[!is.na(ii)]

xx = map$ChEBI_Name %in% ii
tmp = map[xx, ]
tmp = tmp %>% left_join(chebi.charge[, c("ChEBI_ID", "ChEBI_Name", "ChEBI_Name_Convert")], 
                        by = c("ChEBI_Name" = "ChEBI_Name_Convert"))
colnames(tmp) = c("InChikey", "Name", "HMDB_Name", "LM_Name", "ChEBI_Name", "KEGG_Name", 
                  "HMDB", "ChEBI", "KEGG", "LM_ID", "PubChem", "ChEBI_ID_Charge", "ChEBI_Name_Charge")

map = bind_rows(map[!xx, ], tmp)
InChikey_TOP14 = unlist(sapply(map$InChikey, function(x)strsplit(x, "-")[[1]][1]))
map = cbind(map[,1], InChikey_TOP14, map[,2:ncol(map)])
write.csv(map, "compoundID_map.csv", row.names = F, na = "")





















