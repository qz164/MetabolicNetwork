rm(list = ls()[])

library(dplyr)
library(tidyr)

kegg_reactions = read.csv("outputs/kegg_reactions_combined.csv", as.is = T, check.names = F, header = T)
length(unique(kegg_reactions$`Reaction ID`))

##############################################  Update gene symbols for enzyme: Use two annotation files
####################### Uniprot to Gene symbol + Ensembl mapper (BioMart doesn't cover all UniProt)
UniprotHuman = read.delim("inputs/Uniprotkb_Human.tsv", as.is = T, check.names = F, header = T)
UniprotMouse = read.delim("inputs/Uniprotkb_Mouse.tsv", as.is = T, check.names = F, header = T)
UniprotRat = read.delim("inputs/Uniprotkb_Rat.tsv", as.is = T, check.names = F, header = T)

UniprotHuman$Ensembl = gsub(";$", "", UniprotHuman$Ensembl)
UniprotMouse$Ensembl = gsub(";$", "", UniprotMouse$Ensembl)
UniprotRat$Ensembl = gsub(";$", "", UniprotRat$Ensembl)

UniprotHuman$`Gene Names` = unlist(sapply(UniprotHuman$`Gene Names`, function(x)strsplit(x, " ")[[1]][1]))
UniprotMouse$`Gene Names` = unlist(sapply(UniprotMouse$`Gene Names`, function(x)strsplit(x, " ")[[1]][1]))
UniprotRat$`Gene Names` = unlist(sapply(UniprotRat$`Gene Names`, function(x)strsplit(x, " ")[[1]][1]))

UniprotHuman = UniprotHuman %>%
  separate_rows(Ensembl, sep = ";") %>%
  mutate(Ensembl = gsub("\\[.*?\\]", "", Ensembl)) %>%
  mutate(Ensembl = trimws(Ensembl))

UniprotMouse = UniprotMouse %>%
  separate_rows(Ensembl, sep = ";") %>%
  mutate(Ensembl = gsub("\\[.*?\\]", "", Ensembl)) %>%
  mutate(Ensembl = trimws(Ensembl))

UniprotRat = UniprotRat %>%
  separate_rows(Ensembl, sep = ";") %>%
  mutate(Ensembl = gsub("\\[.*?\\]", "", Ensembl)) %>%
  mutate(Ensembl = trimws(Ensembl))

colnames(UniprotHuman)[1] = colnames(UniprotMouse)[1] = colnames(UniprotRat)[1] = "UniProtKB/Swiss-Prot ID"
colnames(UniprotHuman)[9] = colnames(UniprotMouse)[9] = colnames(UniprotRat)[9] = "Transcript stable ID version"

# map by biomart?
mart.human = read.delim("inputs/mart_export_human.txt", as.is = T, check.names = F, header = T)
mart.mouse = read.delim("inputs/mart_export_mouse.txt", as.is = T, check.names = F, header = T)
mart.rat = read.delim("inputs/mart_export_rat.txt", as.is = T, check.names = F, header = T)

MapHuman = UniprotHuman[,c(1,4,9)] %>% 
  full_join(mart.human, by = c("Transcript stable ID version")) %>%
  group_by(`Gene stable ID`) %>%
  summarise(
    `Gene name` = `Gene name`[1],
    `Gene description` = `Gene description`[1],
    Uniprot = paste(na.omit(unique(c(`UniProtKB/Swiss-Prot ID.x`, `UniProtKB/Swiss-Prot ID.y`))), collapse = ",")
  ) %>% 
  separate_rows(Uniprot, sep = ",") %>%
  mutate(
    `Gene name` = ifelse(is.na(`Gene name`), 
                         UniprotHuman$`Gene Names`[match(Uniprot, UniprotHuman$`UniProtKB/Swiss-Prot ID`)],
                         `Gene name`))

MapHuman = MapHuman %>% 
  group_by(`Gene name`) %>%
  summarise(
    `Gene name` = first(na.omit(`Gene name`)),
    `Gene stable ID` = first(na.omit(`Gene stable ID`)),
    `Gene description` = first(na.omit(`Gene description`)),
    Uniprot = paste(unique(na.omit(Uniprot)), collapse = ";")) %>%  
  ungroup() %>% separate_rows(Uniprot, sep = ";") %>%
  group_by(`Gene name`) %>%
  filter(
    !(is.na(Uniprot) | Uniprot == "") |  
      !any(!is.na(Uniprot) & Uniprot != "")) %>% 
  ungroup()


MapMouse = UniprotMouse[,c(1,4,9)] %>% 
  full_join(mart.mouse, by = c("Transcript stable ID version")) %>%
  group_by(`Gene stable ID`) %>%
  summarise(
    `Gene name` = `Gene name`[1],
    `Gene description` = `Gene description`[1],
    Uniprot = paste(na.omit(unique(c(`UniProtKB/Swiss-Prot ID.x`, `UniProtKB/Swiss-Prot ID.y`))), collapse = ",")
  ) %>% 
  separate_rows(Uniprot, sep = ",") %>%
  mutate(
    `Gene name` = ifelse(is.na(`Gene name`), 
                         UniprotMouse$`Gene Names`[match(Uniprot, UniprotMouse$`UniProtKB/Swiss-Prot ID`)],
                         `Gene name`))

MapMouse = MapMouse %>% 
  group_by(`Gene name`) %>%
  summarise(
    `Gene name` = first(na.omit(`Gene name`)),
    `Gene stable ID` = first(na.omit(`Gene stable ID`)),
    `Gene description` = first(na.omit(`Gene description`)),
    Uniprot = paste(unique(na.omit(Uniprot)), collapse = ";")) %>%  
  ungroup() %>% separate_rows(Uniprot, sep = ";") %>%
  group_by(`Gene name`) %>%
  filter(
    !(is.na(Uniprot) | Uniprot == "") |  
      !any(!is.na(Uniprot) & Uniprot != "")) %>% 
  ungroup()


MapRat = UniprotRat[,c(1,4,9)] %>% 
  full_join(mart.rat, by = c("Transcript stable ID version")) %>%
  group_by(`Gene stable ID`) %>%
  summarise(
    `Gene name` = `Gene name`[1],
    `Gene description` = `Gene description`[1],
    Uniprot = paste(na.omit(unique(c(`UniProtKB/Swiss-Prot ID.x`, `UniProtKB/Swiss-Prot ID.y`))), collapse = ",")
  ) %>% 
  separate_rows(Uniprot, sep = ",") %>%
  mutate(
    `Gene name` = ifelse(is.na(`Gene name`), 
                         UniprotRat$`Gene Names`[match(Uniprot, UniprotRat$`UniProtKB/Swiss-Prot ID`)],
                         `Gene name`))

MapRat = MapRat %>% 
  group_by(`Gene name`) %>%
  summarise(
    `Gene name` = first(na.omit(`Gene name`)),
    `Gene stable ID` = first(na.omit(`Gene stable ID`)),
    `Gene description` = first(na.omit(`Gene description`)),
    Uniprot = paste(unique(na.omit(Uniprot)), collapse = ";")) %>%  
  ungroup() %>% separate_rows(Uniprot, sep = ";") %>%
  group_by(`Gene name`) %>%
  filter(
    !(is.na(Uniprot) | Uniprot == "") |  
      !any(!is.na(Uniprot) & Uniprot != "")) %>% 
  ungroup()


####################### Table 1 --> Uniprot mapping From Expasy
enzyme_expasy = read.csv("inputs/enzyme_uniprot_mapping.csv", as.is = T, check.names = F, header = T)

ii = intersect(kegg_reactions$`Participant ID`, enzyme_expasy$`EC Number`)
enzyme_expasy = enzyme_expasy[enzyme_expasy$`EC Number` %in% ii &
                                enzyme_expasy$Organism %in% c("HUMAN", "MOUSE", "RAT"), ]
colnames(enzyme_expasy)[1] = "Participant ID"

# Convert Uniprot ID to gene symbol and ensembl ID
e.human = enzyme_expasy[enzyme_expasy$Organism== "HUMAN",]
e.mouse = enzyme_expasy[enzyme_expasy$Organism == "MOUSE",]
e.rat = enzyme_expasy[enzyme_expasy$Organism == "RAT",]
colnames(MapHuman)[4] = colnames(MapMouse)[4] = colnames(MapRat)[4] = "UniProt ID"

e.human = e.human %>% left_join(MapHuman, by = "UniProt ID")
e.mouse = e.mouse %>% left_join(MapMouse, by = "UniProt ID")
e.rat = e.rat %>% left_join(MapRat, by = "UniProt ID")

enzyme_expasy = rbind(e.human, e.mouse, e.rat)
kegg_reactions = kegg_reactions %>% left_join(enzyme_expasy, by = "Participant ID")

####################### Table 2 --> Gene symbol from KEGG API 
enzyme = read.csv("outputs/enzyme_gene_table.csv", as.is = T, check.names = F, header = T)
colnames(enzyme)[1] = "Participant ID"
enzyme$'Gene name' = gsub("\\)", "", sapply(strsplit(enzyme$gene_symbol, "\\("), `[`, 2))
enzyme$Entrez = sapply(strsplit(enzyme$gene_symbol, "\\("), `[`, 1)
enzyme$Organism = toupper(enzyme$species)
enzyme = enzyme[,-c(2:3)]

ec = unique(kegg_reactions$`Participant ID`[grepl(".",kegg_reactions$`Participant ID`, fixed = T)])
ec.new = NULL
for (i in 1:length(ec)){
  xx = kegg_reactions$`Gene name`[kegg_reactions$`Participant ID` == ec[i]]
  yy = enzyme$`Gene name`[enzyme$`Participant ID` == ec[i]]
  
  id = yy[!yy %in% xx]
  if (length(id) > 0){
    tmp = enzyme[enzyme$`Gene name` %in% id & enzyme$`Participant ID` == ec[i],]
    ec.new = rbind(ec.new, tmp)
  }
}

ec.new = ec.new[!is.na(ec.new$`Gene name`) & ec.new$`Participant ID` %in% kegg_reactions$`Participant ID`,]
df = kegg_reactions[,c("Reaction ID", "Reaction Name", 
                       "Participant ID", "Participant Name", "Type")]
df = df[!duplicated(df),]
ec.new = df %>% left_join(ec.new[,-3], by = c("Participant ID"))
ec.new = ec.new[!is.na(ec.new$Organism),]

####################### Fill the NA cells for genes added from KEGG API
tmp = kegg_reactions %>% full_join(ec.new)
xx = !is.na(tmp$`Gene name`) & is.na(tmp$`Gene stable ID`) & is.na(tmp$`UniProt ID`) & tmp$Organism == "HUMAN"
yy = !is.na(tmp$`Gene name`) & is.na(tmp$`Gene stable ID`) & is.na(tmp$`UniProt ID`) & tmp$Organism == "MOUSE"
zz = !is.na(tmp$`Gene name`) & is.na(tmp$`Gene stable ID`) & is.na(tmp$`UniProt ID`) & tmp$Organism == "RAT"


merged_df = tmp[xx,] %>% left_join(MapHuman, by = c("Gene name")) %>%
  mutate(
    `UniProt ID` = coalesce(`UniProt ID.x`, `UniProt ID.y`),
    `Gene stable ID` = coalesce(`Gene stable ID.x`, `Gene stable ID.y`),
    `Gene description` = coalesce(`Gene description.x`, `Gene description.y`)
  ) %>%
  select(
    `Reaction ID`, `Reaction Name`, `Participant ID`, `Participant Name`, `Type`, 
    `UniProt ID`,`Organism`, `Gene name`, 
    `Gene stable ID`, `Gene description`
  )
tmp[xx,] = merged_df


merged_df = tmp[yy,] %>% left_join(MapMouse, by = c("Gene name")) %>%
  mutate(
    `UniProt ID` = coalesce(`UniProt ID.x`, `UniProt ID.y`),
    `Gene stable ID` = coalesce(`Gene stable ID.x`, `Gene stable ID.y`),
    `Gene description` = coalesce(`Gene description.x`, `Gene description.y`)
  ) %>%
  select(
    `Reaction ID`, `Reaction Name`, `Participant ID`, `Participant Name`, `Type`, 
    `UniProt ID`,`Organism`, `Gene name`, 
    `Gene stable ID`, `Gene description`
  )
tmp[yy,] = merged_df


merged_df = tmp[zz,] %>% left_join(MapMouse, by = c("Gene name")) %>%
  mutate(
    `UniProt ID` = coalesce(`UniProt ID.x`, `UniProt ID.y`),
    `Gene stable ID` = coalesce(`Gene stable ID.x`, `Gene stable ID.y`),
    `Gene description` = coalesce(`Gene description.x`, `Gene description.y`)
  ) %>%
  select(
    `Reaction ID`, `Reaction Name`, `Participant ID`, `Participant Name`, `Type`, 
    `UniProt ID`,`Organism`, `Gene name`, 
    `Gene stable ID`, `Gene description`
  )
tmp[zz,] = merged_df

kegg_reactions = tmp
x = kegg_reactions$`Participant Name` == "N/A"
kegg_reactions$`Participant Name`[x] = kegg_reactions$`Gene name`[x]
kegg_reactions = kegg_reactions[order(kegg_reactions$`Reaction ID`),]


####################### Update compound names
compound = read.delim("inputs/compound.txt", as.is = T, check.names = F, header = F)
colnames(compound) = c("Participant ID", "CompoundName")

kegg_reactions = kegg_reactions %>% left_join(compound, by = "Participant ID")
tmp = unlist(sapply(kegg_reactions$CompoundName, function(x)strsplit(x, ";")[[1]][1]))

x = !is.na(kegg_reactions$CompoundName)
kegg_reactions$`Participant Name`[x] = tmp[x]


####################### Task 1: Any pathway information for each reaction?
pw2reaction = read.delim("inputs/rn.txt", as.is = T, check.names = F, header = F)
pw = read.delim("inputs/pathway.txt", as.is = T, check.names = F, header = F)
colnames(pw) = c("pw", "Name")

colnames(pw2reaction) = c("rn","pw")
pw2reaction = pw2reaction[grep("map", pw2reaction$pw),]
pw2reaction$rn = gsub("rn:", "", pw2reaction$rn)
pw2reaction$pw = gsub("path:", "", pw2reaction$pw)

m = match(pw2reaction$pw, pw$pw)
pw2reaction$PathwayName = pw$Name[m]

m = match(kegg_reactions$`Reaction ID`, pw2reaction$rn)
kegg_reactions = cbind("PathwayID" = pw2reaction$pw[m], "Pathway" = pw2reaction$PathwayName[m],
                       kegg_reactions)

# update reaction name
reaction = read.delim("inputs/reaction.txt", as.is = T, check.names = F, header = F)
colnames(reaction) = c("rn", "Name")

m = match(kegg_reactions$`Reaction ID`, reaction$rn)
kegg_reactions$`Reaction Name` = reaction$Name[m]
#write.table(kegg_reactions, "kegg_reactions_updated.txt", na = "", quote = F, sep = "\t", row.names = F)


####################### Task2: Identify common metabolites
cpd.reactions = kegg_reactions[grepl("C", kegg_reactions$`Participant ID`) &
                                 ! grepl(".", kegg_reactions$`Participant ID`, fixed = T), ]

tmp = cpd.reactions[,c("Participant ID", "Participant Name")]
tmp = tmp[!duplicated(tmp),]

tmp$Freq = NA
for (i in 1:nrow(tmp)){
  tmp$Freq[i] = length(unique(cpd.reactions$`Reaction ID`[cpd.reactions$`Participant ID` == tmp$`Participant ID`[i]]))
}

m = match(kegg_reactions$`Participant ID`, tmp$`Participant ID`)
kegg_reactions$CompoundFreq = tmp$Freq[m]


####################### Task3: Add mapping ID
kegg_reactions$KEGG_ID = ifelse(grepl("C", kegg_reactions$`Participant ID`) & 
                                    !grepl(".", kegg_reactions$`Participant ID`, fixed = T),
                                  kegg_reactions$`Participant ID`, NA)

####################### Task4: separate human and mouse?
library(dplyr)

kegg.human = kegg_reactions %>%
  group_by(`Reaction ID`) %>%
  filter(any(grepl("HUMAN", Organism))) %>%
  filter(!grepl("MOUSE|RAT", Organism)) %>%
  filter(!(is.na(`Participant Name`))) %>%
  ungroup()

kegg.mouse = kegg_reactions %>%
  group_by(`Reaction ID`) %>%
  filter(any(grepl("MOUSE", Organism))) %>%
  filter(!grepl("HUMAN|RAT", Organism)) %>%
  filter(!(is.na(`Participant Name`))) %>%
  ungroup()

kegg.rat = kegg_reactions %>%
  group_by(`Reaction ID`) %>%
  filter(any(grepl("RAT", Organism))) %>%
  filter(!grepl("HUMAN|MOUSE", Organism)) %>%
  filter(!(is.na(`Participant Name`))) %>%
  ungroup()

write.table(kegg.human, "final_parsed_info/KEGG_Reactions_HSA.txt", sep = "\t", quote = F, na = "", row.names = F)
write.table(kegg.mouse, "final_parsed_info/KEGG_Reactions_MMU.txt", sep = "\t", quote = F, na = "", row.names = F)
write.table(kegg.rat, "final_parsed_info/KEGG_Reactions_RNO.txt", sep = "\t", quote = F, na = "", row.names = F)

# Count metabolites and genes
tmp = kegg.human
length(unique(tmp$KEGG_ID))
length(unique(tmp$`Gene stable ID`))





