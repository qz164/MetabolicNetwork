# R/build_reactome.R

#' Constructs a metabolic network from Reactome data for specified species.
#'
#' @param Merged_data Input dataset containing merged metabolic and transcriptomic data
#' @param info Metadata for column standardization
#' @param species Species identifier ("rno", "mmu", or "hsa")
#' @param gene_id Column name for gene identifiers (default: "Uniprot")
#' @param met_id Column name for metabolite identifiers (default: "ChEBI")
#' @param size Network size type (default: "MetMN_Expand")
#' @param write_output Should output files be written? (default: FALSE)
#' @param output_dir Output directory path (default: "outputs")
#' @return List containing node attributes, edge attributes, and solid edges
#' @export
#' @import dplyr
#' @import tidyr
#' @import purrr
#' @import stringr
#' @importFrom utils write.table
#'
build_reactome <- function(Merged_data,
                           info,
                           species = "rno",
                           gene_id,
                           met_id,
                           size = "MetMN_Expand",
                           write_output = FALSE,
                           output_dir = "outputs") {

  # Load species-specific data (INTERNAL)
  reaction <- switch(
    species,
    "rno" = MetabolicNetwork:::reactome.reaction.rno,
    "mmu" = MetabolicNetwork:::reactome.reaction.mmu,
    "hsa" = MetabolicNetwork:::reactome.reaction.hsa,
    stop("Invalid species: ", species)
  )

  # Use other internal datasets directly
  ont_chebi <- MetabolicNetwork:::ont_chebi
  chebi_data <- MetabolicNetwork:::chebi
  map <- MetabolicNetwork:::map

  # Core processing
  result <- .build_reactome_core(
    Merged_data = Merged_data,
    info = info,
    reaction = reaction,
    ont_chebi = ont_chebi,
    chebi_data = chebi,
    map = map,
    gene_id = gene_id,
    met_id = met_id,
    size = size
  )

  # Write outputs
  if (write_output) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  utils::write.table(result$nodeAttrib,
                    file.path(output_dir, "Reactome_nodes.txt"),
                    sep = "\t", quote = FALSE,
                    na = "", row.names = FALSE)

  utils::write.table(result$edgeAttrib,
                    file.path(output_dir, "Reactome_edges.txt"),
                    sep = "\t", quote = FALSE,
                    na = "", row.names = FALSE)

  utils::write.table(result$edgeAttrib_solid,
                    file.path(output_dir, "Reactome_edges_solid.txt"),
                    sep = "\t", quote = FALSE,
                    na = "", row.names = FALSE)
  }
}


.build_reactome_core <- function(Merged_data,
                                 info,
                                 reaction,
                                 ont_chebi,
                                 chebi_data,
                                 map,
                                 gene_id,
                                 met_id,
                                 size) {

  ################## Cleanup identifies
  exp.data = Merged_data
  id = ifelse(!is.na(exp.data[[gene_id]]), exp.data[[gene_id]],
              ifelse(!is.na(exp.data[[met_id]]), exp.data[[met_id]], NA))
  exp.data = cbind("MappingID" = id, exp.data)
  exp.data = exp.data[!is.na(exp.data$MappingID),]

  # Remove duplications: for the same mapping ID,
  # choose the most representative compound/genes (maximal absolute significance)
  exp.data = exp.data %>%
    group_by(MappingID) %>%
    filter(if (any(`Differential feature` == '*' & !is.na(`Differential feature`))) {
      majority_sign = ifelse(sum(sign(MaxSig) > 0) >= sum(sign(MaxSig) < 0), 1, -1)
      `Differential feature` == '*' & sign(MaxSig) == majority_sign & abs(MaxSig) == max(abs(MaxSig))
    } else {
      abs(MaxSig) == max(abs(MaxSig))
    }) %>%
    ungroup() %>%
    group_by(MappingID) %>%
    dplyr::slice(1) %>%
    ungroup()

  ################## Clean up Reactome data
  ## ID mapping: the reaction must contain at least one measured molecule
  reaction[is.na(reaction)] = ""
  reaction$ChEBI_ID = as.character(reaction$ChEBI_ID)

  ## Remove monoatomic ion in the database
  r = c("is_a", "is_conjugate_base_of",
        "is_enantiomer_of", "has_functional_parent", "is_conjugate_acid_of",
        "is_tautomer_of")
  ont_chebi = ont_chebi[ont_chebi$Relationship %in% r, ]


  ion = .find_offspring("CHEBI:24867",
                       ont_chebi, visited = NULL) # remove all ions (monoatomic ion)
  rid = chebi$ChEBI_Name[chebi$ChEBI_ID %in% unlist(ion)] # check names
  ion = gsub("CHEBI:", "", unlist(ion))
  reaction = reaction[!reaction$ChEBI_ID %in% ion, ]

  ## Remove metabolites with low weight (<50)
  mw = chebi
  mw$ChEBI_ID = gsub("CHEBI:", "", mw$ChEBI_ID)
  low_mw = mw$ChEBI_ID[mw$`Monoisotopic Mass` < 50]
  low_mw = low_mw[low_mw != ""]
  reaction = reaction[!reaction$ChEBI_ID %in% low_mw , ]


  ######################################### Network construction ########################################################
  ### Step 1: Remove reactions with no entities detected. Choose "Metabolism" only?
  network1 = reaction %>% group_by(`Reaction ID`) %>%
    filter(any(ChEBI_ID %in% exp.data$MappingID) |
             any(UniProt %in% exp.data$MappingID))  %>%
    filter(grepl("Metabolism", TopPathway)) # Metabolic Pathways only

  # Label those metabolites that were detected
  network1$ChEBI_ID_exp = ifelse(network1$ChEBI_ID %in% exp.data$MappingID, "*", "")
  network1$UniProt_exp = ifelse(network1$UniProt %in% exp.data$MappingID, "*", "")

  ### Step2: Retrieve and filter pair relationships from reactions
  nodes = exp.data$MappingID[exp.data$`Differential feature` == "*" & !is.na(exp.data$`Differential feature`)]
  nodes = unique(nodes[nodes != ""])

  tmp = network1 %>% group_by(`Reaction ID`) %>%
    mutate(SeedNodes = ifelse(ChEBI_ID %in% nodes | UniProt %in% nodes, "*", "")) %>%
    mutate(MappingID = ifelse(ChEBI_ID == "", UniProt, ChEBI_ID)) %>%
    mutate(Detected = ifelse(ChEBI_ID_exp == "*" | UniProt_exp == "*", "*", "")) %>%
    filter(any(ChEBI_ID %in% nodes) | any(UniProt %in% nodes)) %>%
    filter(n() > 1) %>% ungroup()

  m = match(tmp$MappingID, exp.data$MappingID)
  tmp$`Differential feature` = exp.data$`Differential feature`[m]
  tmp[is.na(tmp)] = ""
  tmp$molecule = ifelse(grepl("[A-Za-z]", tmp$MappingID), "GeneProduct", "Metabolite")


  # Transform into pairs
  pairs_df = tmp %>%
    group_by(`Reaction ID`) %>%
    summarise(
      pairs = list(combn(participant_name, 2, simplify = FALSE)),
      types = list(combn(Type, 2, simplify = FALSE)),
      seeds = list(combn(SeedNodes, 2, simplify = FALSE)),
      mappingID = list(combn(MappingID, 2, simplify = FALSE)),
      detectedPairs = list(combn(Detected, 2, simplify = FALSE)),
      Entity = list(combn(`Schema Class`, 2, simplify = FALSE)),
      Relationship = list(combn(molecule, 2, simplify = FALSE)),
      Sig = list(combn(`Differential feature`, 2, simplify = FALSE)),
      .groups = 'drop'
    ) %>%
    unnest(c(pairs, types, seeds, mappingID, detectedPairs, Sig, Entity, Relationship)) %>%
    mutate(pair1 = map_chr(pairs, 1),
           pair2 = map_chr(pairs, 2),
           type1 = map_chr(types, 1),
           type2 = map_chr(types, 2),
           SeedNode1 = map_chr(seeds, 1),
           SeedNode2 = map_chr(seeds, 2),
           mappingID1 = map_chr(mappingID, 1),
           mappingID2 = map_chr(mappingID, 2),
           detectedPair1 = map_chr(detectedPairs, 1),
           detectedPair2 = map_chr(detectedPairs, 2),
           SigNode1 = map_chr(Sig, 1),
           SigNode2 = map_chr(Sig, 2),
           Entity1 = map_chr(Entity, 1),
           Entity2 = map_chr(Entity, 2),
           Relationship = paste0(map_chr(Relationship, 1), "_", map_chr(Relationship, 2)),
           #Timepoint1 = map_chr(Timepoint, 1),
           #Timepoint2 = map_chr(Timepoint, 2),
    ) %>%
    filter(mappingID1 != mappingID2 & (SeedNode1 == "*" | SeedNode2 == "*")) %>%
    filter(mappingID1 != "" & mappingID2 != "") %>%
    ungroup() %>%
    mutate(sorted_pairs = paste0(pmin(mappingID1, mappingID2), "_", pmax(mappingID1, mappingID2))) %>%
    mutate(Relationship = gsub("GeneProduct_Metabolite", "Metabolite_GeneProduct", Relationship)) %>%
    dplyr::select(-pairs, -types, -seeds, -mappingID, -detectedPairs, -Sig, -Entity)

  complex = c("CandidateSet", "Complex", "DefinedSet")
  pairs_df = pairs_df %>%
    filter(!(type1 == "catalyst" & type2 == "catalyst")) %>%
    filter(!(Relationship == "GeneProduct_GeneProduct" & type1 != type2 & all(c(type1, type2) != "catalyst") &
               Entity1 %in% complex & Entity2 %in% complex)) %>% # Remove --> two proteins from input-complex and output-complex
    mutate(EdgeColor = ifelse(type1 == type2, "SameComplex", "")) %>%
    distinct()

  # Apply aggregation to concatenate unique values for each group
  network2 = aggregate(. ~ type1 + type2 + mappingID1 + mappingID2,
                       data = pairs_df, FUN = .concat_unique)


  # For non-Sig nodes, only keep it if it connect at least two sig nodes
  nonSig = unique(tmp$MappingID[tmp$`Differential feature` != "*"])
  rid = NULL
  for (i in 1:length(nonSig)){
    node = nonSig[i]
    tmp = network2[network2$mappingID1 == node | network2$mappingID2 == node, ]
    if (nrow(tmp) == 1){
      rid = c(rid, node)
    }
  }
  network2 = network2[!network2$mappingID1 %in% rid & !network2$mappingID2 %in% rid, ]

  ### Step3. keep core network and nodes linked to it

  # The effect sizes can be very different across different studies
  # We thus provide 3 options to visualize networks at different resolution
  # * means differential features

  ## Metabolite-centric metabolic network for differential features only ("MetMN_Diff")
  # For studies with large numbers of differential features
  # Nodes are all differential features and with linkage to metabolites

  ## Metabolite-centric metabolic network with expansion ("MetMN_Expand")
  # For studies with less differential feature
  # an investigation of their neighbors in the network
  # Diff metabolites as seeds ---> *met-*met and *met-*gene
  # First neighbors of the seeds; preserve relationships between first neighbors as well

  ## Metabolic network for differential features ("MN_Diff")
  # To visualize the relationships between all differential features
  # (both metabolites and genes)
  network3 = .clean_metabolic_network_reactome(network2, size)
  print("Network built")

  ###############################
  ###### Prepare node attribute file
  ###############################
  nodeAttrib = data.frame(
    MappingID = c(network3$mappingID1, network3$mappingID2),
    Label_db = c(network3$pair1, network3$pair2),
    detected = c(network3$detectedPair1, network3$detectedPair2),
    SigNode = c(network3$SigNode1, network3$SigNode2)
  )
  nodeAttrib = nodeAttrib[!duplicated(nodeAttrib),]
  nodeAttrib$NodeType = ifelse(grepl("[A-Za-z]", nodeAttrib$MappingID),
                               "GeneProduct", "Metabolite")

  c = info$Colname
  m = match(nodeAttrib$MappingID, exp.data$MappingID)
  nodeAttrib = cbind(nodeAttrib, exp.data[m,c])

  nodeAttrib$Label_db = unlist(lapply(nodeAttrib$Label_db, function(x)strsplit(x, " [", fixed = T)[[1]][1]))
  x = nodeAttrib$detected == ""
  nodeAttrib$Name_map[x] = nodeAttrib$Label_db[x] # For undetected genes --> (Only uniprot and gene symbol available)

  ####### Meta data for undetected nodes
  # For undetected metabolites --> (ChEBI available; get HMDB ID, KEGG)
  x = nodeAttrib$NodeType == "Metabolite" & nodeAttrib$detected == ""
  nodeAttrib.unmapp = nodeAttrib[x,]

  map$ChEBI = as.character(map$ChEBI)
  map$ChEBI_ID_Charge = as.character(map$ChEBI_ID_Charge)

  nodeAttrib.unmapp = nodeAttrib.unmapp %>%
    left_join(map[,c("ChEBI", "KEGG", "HMDB")], by = c("MappingID" = "ChEBI")) %>%
    mutate(
      KEGG = coalesce(KEGG.x, KEGG.y),
      `HMDB lip` = coalesce(HMDB, `HMDB lip`),
      `HMDB accession` = coalesce(HMDB, `HMDB accession`)
    ) %>%
    dplyr::select(-KEGG.x, -KEGG.y, -HMDB)

  nodeAttrib.unmapp2 = nodeAttrib.unmapp %>% filter(is.na(KEGG) & is.na(`HMDB lip`))
  nodeAttrib.unmapp = nodeAttrib.unmapp %>% filter(!is.na(KEGG) | !is.na(`HMDB lip`))

  # Second match by ChEBI_ID_Charge only for unmatched
  nodeAttrib.unmapp2 = nodeAttrib.unmapp2 %>%
    left_join(map[, c("ChEBI_ID_Charge", "KEGG", "HMDB")], by = c("MappingID" = "ChEBI_ID_Charge")) %>%
    mutate(
      KEGG = coalesce(KEGG.x, KEGG.y),
      `HMDB lip` = coalesce(HMDB, `HMDB lip`),
      `HMDB accession` = coalesce(HMDB, `HMDB accession`)
    ) %>%
    select(-KEGG.x, -KEGG.y,  -HMDB)
  nodeAttrib = bind_rows(nodeAttrib[!x,], nodeAttrib.unmapp, nodeAttrib.unmapp2)
  # This didn't solve the problem for lipid names

  # Compartment
  tmp = pairs_df %>%
    filter(sorted_pairs %in% network3$sorted_pairs)
  comp = data.frame("Reaction ID" = c(tmp$`Reaction ID`, tmp$`Reaction ID`),
                    "Compartment" = c(tmp$pair1, tmp$pair2),
                    "MappingID" = c(tmp$mappingID1, tmp$mappingID2))

  comp$Compartment = sapply(comp$Compartment,
                            function(x)trimws(str_extract_all(x, " \\[([^\\]]+)\\]")[[1]]))
  comp = comp[!duplicated(comp),]
  nodeAttrib = nodeAttrib %>% left_join(comp, by = "MappingID")
  nodeAttrib[is.na(nodeAttrib)] = ""

  # label significant nodes only
  diff = nodeAttrib[,grepl("_diff", colnames(nodeAttrib))]
  seq = info$Type[match(colnames(diff), info$Colname)]
  seq = as.integer(gsub("Condition ", "", seq))

  tmp = matrix(ncol = length(seq), nrow = nrow(nodeAttrib), "")
  colnames(tmp) = paste0("Label_Sig", seq)

  for (i in 1:ncol(tmp)){
    tmp[,i] = ifelse(diff[,i] == "*", nodeAttrib$Name_map, "")
  }
  nodeAttrib = cbind(nodeAttrib, tmp)
  nodeAttrib$Label_Sig = ifelse(nodeAttrib$`Differential feature` == "*", nodeAttrib$Name_map, "")
  nodeAttrib = aggregate(. ~ MappingID,
                         data = nodeAttrib, FUN = .concat_unique)

  ############ Node annotation
  # Log2 Fold change as color
  col = info[!is.na(info$NodeColor),]
  m = match(col$Colname, colnames(exp.data))
  tmp = exp.data[,m]
  colnames(tmp) = col$NodeColor

  m = match(nodeAttrib$MappingID, exp.data$MappingID)
  tmp = tmp[m, ]
  tmp[is.na(tmp)] = 0
  nodeAttrib = cbind(nodeAttrib, tmp)

  # Max absolute value
  nodeAttrib$ColorAny = apply(tmp, 1, function(x)x[which.max(abs(x))])


  # Significance as size
  size = info[!is.na(info$NodeSize),]
  m = match(size$Colname, colnames(exp.data))
  tmp = exp.data[,m]
  colnames(tmp) = size$NodeSize

  m = match(nodeAttrib$MappingID, exp.data$MappingID)
  tmp = tmp[m, ]
  tmp[is.na(tmp)] = 0.1
  tmp = abs(tmp)
  tmp[tmp <= 0.1] = 0.1
  nodeAttrib = cbind(nodeAttrib, tmp)

  # Max absolute value
  nodeAttrib$SizeAny = apply(tmp, 1, function(x)x[which.max(x)])
  print("Node attribute file")

  ###############################
  ###### Prepare edge attribute file
  ###############################
  network3$Direction = paste0(network3$type1, "_", network3$type2)

  network3 = network3 %>%
    group_by(sorted_pairs) %>%
    filter(if (all(EdgeColor == "SameComplex")) row_number() == 1 else TRUE) %>%
    ungroup()

  # Highlight those connections with concordant changes between enzyme and substrate/product
  concordant = rep("", nrow(network3))
  for (i in 1:nrow(network3)){
    r = as.character(network3[i, c("type1", "type2")])
    catalyst = r == "catalyst"
    prot = as.character(network3[i, c("mappingID1", "mappingID2")[catalyst]])
    met = as.character(network3[i, c("mappingID1", "mappingID2")[!catalyst]])

    if (any(catalyst) & length(prot) > 0 & length(met) > 0){
      if (grepl("[A-Za-z]",prot) &
          !grepl("[A-Za-z]", met)){

        if ("output" %in% r){
          m = match(c(prot, met), nodeAttrib$MappingID)
          tmp = as.numeric(nodeAttrib$MaxSig[m])
          xx = (all(sign(tmp) == 1) | all(sign(tmp) == -1)) & all(network3[i, c("SigNode1", "SigNode2")] == "*")
          concordant[i] = ifelse(xx, "Enzyme-Product",concordant[i])
        }

        if ("input" %in% r){
          m = match(c(prot, met), nodeAttrib$MappingID)
          tmp = as.numeric(nodeAttrib$MaxSig[m])
          xx = all(c(-1, 1) %in% sign(tmp)) & all(network3[i, c("SigNode1", "SigNode2")] == "*")
          concordant[i] = ifelse(xx, "Substrate-Enzyme",concordant[i])
        }

        network3$EdgeColor[i] = "Enzyme or transporter"
      }
    }
  }
  network3 = cbind(network3, concordant)
  network3$EdgeColor = ifelse(network3$EdgeColor == "" & network3$Relationship == "Metabolite_Metabolite"
                              & !grepl("catalyst", network3$Direction), "Metabolic Transformation", network3$EdgeColor)
  network3$EdgeColor[network3$EdgeColor == ""] = "Others"
  network3$EdgeLine = ifelse(network3$detectedPair1 == "*" & network3$detectedPair2 == "*", "Solid", "Dashed")

  network3 = network3 %>% group_by(sorted_pairs, EdgeColor) %>%
    dplyr::slice(1)

  ####### Remove duplicated edges
  priority_order = c("Enzyme or transporter", "Metabolic Transformation", "SameComplex", "Others")
  network3 = network3 %>%
    group_by(sorted_pairs) %>%
    mutate(
      EdgeColor = if (n() > 1) {
        tmp = as.character(EdgeColor)
        tmp[which.min(match(tmp, priority_order))]
      } else {
        EdgeColor
      }
    ) %>%
    summarise(
      ReactionID = paste(unique(na.omit(na_if(`Reaction ID`, ""))), collapse = ", "),  # Summarize Reaction ID
      across(-ReactionID, first)
    ) %>%
    ungroup()


  ############ Highlight edges annotation connecting two significant nodes
  m1 = match(network3$mappingID1, nodeAttrib$MappingID)
  m2 = match(network3$mappingID2, nodeAttrib$MappingID)

  diff_cols = grep("_diff", names(nodeAttrib), value = TRUE)

  for (col in diff_cols) {
    network3[[col]] = ifelse(
      nodeAttrib[[col]][m1] == "*" | nodeAttrib[[col]][m2] == "*",
      "*", ""
    )
  }
  edgeAttrib = network3
  edgeAttrib_solid = network3 %>% filter(EdgeLine == "Solid")
  print("Edge attribute file")

  return(list(
    nodeAttrib = nodeAttrib,
    edgeAttrib = edgeAttrib,
    edgeAttrib_solid = edgeAttrib_solid
  ))
}











