# R/build_kegg.R

#' Constructs a metabolic network from KEGG data for specified species.
#'
#' @param Merged_data Input dataset containing merged metabolic and transcriptomic data
#' @param info Metadata for column standardization
#' @param species Species identifier ("rno", "mmu", or "hsa")
#' @param gene_id Column name for gene identifiers
#' @param met_id Column name for metabolite identifiers
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
build_kegg <- function(Merged_data,
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
    "rno" = MetabolicNetwork:::kegg.reaction.rno,
    "mmu" = MetabolicNetwork:::kegg.reaction.mmu,
    "hsa" = MetabolicNetwork:::kegg.reaction.hsa,
    stop("Invalid species: ", species)
  )

  # Use other internal datasets directly
  map <- MetabolicNetwork:::map

  # Core processing
  result <- .build_kegg_core(
    Merged_data = Merged_data,
    info = info,
    reaction = reaction,
    map = map,
    gene_id = gene_id,
    met_id = met_id,
    size = size
  )

  # Write outputs
  if (write_output) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

    utils::write.table(result$nodeAttrib,
                       file.path(output_dir, "KEGG_nodes.txt"),
                       sep = "\t", quote = FALSE,
                       na = "", row.names = FALSE)

    utils::write.table(result$edgeAttrib,
                       file.path(output_dir, "KEGG_edges.txt"),
                       sep = "\t", quote = FALSE,
                       na = "", row.names = FALSE)

    utils::write.table(result$edgeAttrib_solid,
                       file.path(output_dir, "KEGG_edges_solid.txt"),
                       sep = "\t", quote = FALSE,
                       na = "", row.names = FALSE)
  }
}
.build_kegg_core <- function(Merged_data,
                                 info,
                                 reaction,
                                 map,
                                 gene_id,
                                 met_id,
                                 size) {
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


  ########################### Metabolic network
  ##### ID mapping: the reaction must contain at least one measured molecule
  reaction[is.na(reaction)] = ""
  exp.data = as.data.frame(exp.data)
  exp.data[is.na(exp.data)] = ""

  # Filter 1: at least one element was detected in our data.
  network1 = reaction %>% group_by(`Reaction ID`) %>%
    filter(any(KEGG_ID %in% exp.data$MappingID) |
             any(`Gene name` %in% exp.data$MappingID))

  # Label those metabolites that were detected
  network1$KEGG_ID_exp = ifelse(network1$KEGG_ID %in% exp.data$MappingID, "*", "")
  network1$Gene_exp = ifelse(network1$`Gene name` %in% exp.data$MappingID, "*", "")

  # seed nodes: significant metabolites or significant proteins
  nodes = exp.data$MappingID[exp.data$`Differential feature` == "*" & !is.na(exp.data$`Differential feature`)]
  nodes = unique(nodes[nodes != ""])

  # Remove common metabolites with very low weight
  network1 = network1[!network1$`Participant ID` %in% c("C00001", "C00080", "C00007", "C00011"), ]

  tmp = network1 %>% group_by(`Reaction ID`) %>%
    mutate(SeedNodes = ifelse(KEGG_ID %in% nodes | `Gene name` %in% nodes, "*", "")) %>%
    mutate(MappingID = ifelse(KEGG_ID == "", `Gene name`, KEGG_ID)) %>%
    mutate(Detected = ifelse(KEGG_ID_exp == "*" | Gene_exp == "*", "*", "")) %>%
    filter(any(KEGG_ID %in% nodes) | any(`Gene name` %in% nodes)) %>%
    filter(n() > 1) %>% ungroup()

  m = match(tmp$MappingID, exp.data$MappingID)
  tmp$`Differential feature` = exp.data$`Differential feature`[m]
  #tmp$FirstResponse = as.character(exp.data$FirstResponse[m])
  tmp[is.na(tmp)] = ""
  tmp$molecule = ifelse(grepl("^C\\d{5}$", tmp$MappingID), "Metabolite", "GeneProduct")

  pairs_df = tmp %>%
    group_by(`Reaction ID`) %>%
    summarise(pairs = list(combn(`Participant Name`, 2, simplify = FALSE)),
              types = list(combn(Type, 2, simplify = FALSE)),
              seeds = list(combn(SeedNodes, 2, simplify = FALSE)),
              mappingID = list(combn(MappingID, 2, simplify = FALSE)),
              detectedPairs = list(combn(Detected, 2, simplify = FALSE)),
              Sig = list(combn(`Differential feature`, 2, simplify = FALSE)),
              Relationship = list(combn(molecule, 2, simplify = FALSE)),
              .groups = 'drop') %>%
    unnest(c(pairs, types, seeds, mappingID, detectedPairs, Sig, Relationship)) %>%
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
           Relationship = paste0(map_chr(Relationship, 1), "_", map_chr(Relationship, 2)),
    ) %>%
    filter(type1 != type2 & pair1 != pair2 & mappingID1 != mappingID2 &
             (SeedNode1 == "*" | SeedNode2 == "*")) %>%
    filter(mappingID1 != "" & mappingID2 != "") %>%
    ungroup() %>%
    mutate(sorted_pairs = paste0(pmin(mappingID1, mappingID2), "_", pmax(mappingID1, mappingID2))) %>%
    dplyr::select(-pairs, -types, -seeds, -mappingID, -detectedPairs, -Sig)

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
  network3 = .clean_metabolic_network_kegg(network2, size)
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
  nodeAttrib$NodeType = ifelse(grepl("^C\\d{5}$", nodeAttrib$MappingID),
                               "Metabolite", "GeneProduct")

  c = info$Colname
  m = match(nodeAttrib$MappingID, exp.data$MappingID)
  nodeAttrib = cbind(nodeAttrib, exp.data[m,c])

  x = nodeAttrib$detected == ""
  nodeAttrib$Name_map[x] = nodeAttrib$Label_db[x] # For undetected genes --> (Only uniprot and gene symbol available)

  ####### Meta data for undetected nodes
  x = nodeAttrib$NodeType == "Metabolite" & nodeAttrib$detected == ""
  nodeAttrib.unmapp = nodeAttrib[x,]

  map$ChEBI = as.character(map$ChEBI)
  nodeAttrib.unmapp = nodeAttrib.unmapp %>%
    left_join(map[,c("ChEBI", "KEGG", "HMDB")], by = c("MappingID" = "KEGG")) %>%
    mutate(
      ChEBI = coalesce(ChEBI.x, ChEBI.y),
      `HMDB lip` = coalesce(HMDB, `HMDB lip`),
      `HMDB accession` = coalesce(HMDB, `HMDB accession`)
    ) %>%
    dplyr::select(-ChEBI.x, -ChEBI.y, -HMDB)

  nodeAttrib = bind_rows(nodeAttrib[!x,], nodeAttrib.unmapp)

  tmp = pairs_df %>%
    filter(sorted_pairs %in% network3$sorted_pairs)
  tmp = data.frame("Reaction ID" = c(tmp$`Reaction ID`, tmp$`Reaction ID`),
                   "MappingID" = c(tmp$mappingID1, tmp$mappingID2))
  tmp = tmp[!duplicated(tmp),]
  nodeAttrib = nodeAttrib %>% left_join(tmp, by = "MappingID")

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

  nodeAttrib[is.na(nodeAttrib)] = ""
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
  tmp = apply(tmp, 2, as.numeric)
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
  tmp = apply(tmp, 2, as.numeric)
  tmp = abs(tmp)
  tmp[tmp <= 0.1] = 0.1
  nodeAttrib = cbind(nodeAttrib, tmp)

  # Max absolute value
  nodeAttrib$SizeAny = apply(tmp, 1, function(x)x[which.max(x)])
  print("Node attribute file")

  ###############################
  ###### Prepare edge attribute file
  ###############################
  network3$EdgeColor = ifelse(network3$Relationship == "Metabolite_GeneProduct",
                              "Enzyme", "Metabolic Transformation")
  network3$EdgeLine = ifelse(network3$detectedPair1 == "*" & network3$detectedPair2 == "*", "Solid", "Dashed")

  network3 = network3 %>%  group_by(sorted_pairs) %>%
    filter(row_number() == 1)

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
