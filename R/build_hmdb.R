# R/build_hmdb.R

#' Constructs a metabolic network from HMDB data for specified species.
#'
#' @param Merged_data Input dataset containing merged metabolic and transcriptomic data
#' @param info Metadata for column standardization
#' @param species Species identifier ("rno", "mmu", or "hsa")
#' @param gene_id Column name for gene identifiers
#' @param met_id Column name for metabolite identifiers
#' @param size Network size type (default: "Diff_Expand")
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
build_hmdb <- function(Merged_data,
                       info,
                       species = "rno",
                       gene_id,
                       met_id,
                       size = "Diff_Expand",
                       write_output = FALSE,
                       output_dir = "outputs") {

  # Load species-specific data (INTERNAL)
  HMDB_pair <- switch(
    species,
    "rno" = MetabolicNetwork:::HMDB.reaction.rno,
    "mmu" = MetabolicNetwork:::HMDB.reaction.mmu,
    "hsa" = MetabolicNetwork:::HMDB.reaction.hsa,
    stop("Invalid species: ", species)
  )

  # Use other internal datasets directly
  map <- MetabolicNetwork:::map

  # Core processing
  result <- .build_hmdb_core(
    Merged_data = Merged_data,
    info = info,
    HMDB_pair = HMDB_pair,
    map = map,
    gene_id = gene_id,
    met_id = met_id,
    size = size
  )

  # Write outputs
  if (write_output) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

    utils::write.table(result$nodeAttrib,
                       file.path(output_dir, "HMDB_nodes.txt"),
                       sep = "\t", quote = FALSE,
                       na = "", row.names = FALSE)

    utils::write.table(result$edgeAttrib,
                       file.path(output_dir, "HMDB_edges.txt"),
                       sep = "\t", quote = FALSE,
                       na = "", row.names = FALSE)

    utils::write.table(result$edgeAttrib_solid,
                       file.path(output_dir, "HMDB_edges_solid.txt"),
                       sep = "\t", quote = FALSE,
                       na = "", row.names = FALSE)
  }
}
.build_hmdb_core <- function(Merged_data,
                                 info,
                                 HMDB_pair,
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


  ######################## Metabolic Network
  sig = exp.data$MappingID[exp.data$`Differential feature` == "*"]

  core = HMDB_pair %>%
    filter(metabolite_accession %in% sig & GeneName %in% sig)

  # preserve the linking nodes that connect to those core edges
  sig = unique(c(core$metabolite_accession, core$GeneName))

  if (size == "Diff_Expand"){
    network = HMDB_pair[HMDB_pair$metabolite_accession %in% sig | HMDB_pair$GeneName %in% sig, ]
  } else if (size == "Diff") {
    network = HMDB_pair[HMDB_pair$metabolite_accession %in% sig & HMDB_pair$GeneName %in% sig, ]
  }

  # for those non sig --> must connect with two sig
  nonsig = c(network$metabolite_accession, network$GeneName)
  nonsig = nonsig[!nonsig %in% sig]
  rid = nonsig[!(duplicated(nonsig, fromLast = F) | duplicated(nonsig, fromLast = T))]
  network2 = network[!(network$metabolite_accession %in% rid | network$GeneName %in% rid), ]

  network2$detectedPair1 = ifelse(network2$metabolite_accession %in% exp.data$MappingID, "*", "")
  network2$detectedPair2 = ifelse(network2$GeneName %in% exp.data$MappingID, "*", "")

  network2$EdgeColor = network2$protein_type
  network2$EdgeLine = ifelse(network2$detectedPair1 == "*" & network2$detectedPair2 == "*", "Solid", "Dashed")

  network2 = network2 %>%
    mutate(sorted_pairs = paste0(metabolite_accession, "_", GeneName))

  colnames(network2)[colnames(network2) == "metabolite_accession"] = "mappingID1"
  colnames(network2)[colnames(network2) == "GeneName"] = "mappingID2"
  print("Network built")

  ###############################
  ###### Prepare node attribute file
  ###############################
  nodeAttrib = data.frame(
    MappingID = c(network2$mappingID1, network2$mappingID2),
    name_HMDB = c(network2$MetaboliteName, network2$mappingID2))
  nodeAttrib = nodeAttrib[!duplicated(nodeAttrib),]
  nodeAttrib$NodeType = ifelse(grepl("HMDB", nodeAttrib$MappingID), "Metabolite", "GeneProduct")
  colnames(nodeAttrib)[colnames(nodeAttrib) == "name_HMDB"] = "Label_db"

  c = info$Colname
  m = match(nodeAttrib$MappingID, exp.data$MappingID)
  nodeAttrib = cbind(nodeAttrib, exp.data[m,c])
  #nodeAttrib[is.na(nodeAttrib)] = ""

  nodeAttrib$detected = ifelse(!is.na(nodeAttrib$Name_map), "*", "") # Check
  x = nodeAttrib$detected == ""
  nodeAttrib$Name_map[x] = nodeAttrib$Label_db[x] # For undetected genes --> (Only uniprot and gene symbol available)

  ####### Meta data for undetected nodes
  x = nodeAttrib$NodeType == "Metabolite" & nodeAttrib$detected == ""
  nodeAttrib.unmapp = nodeAttrib[x,]

  #map$ChEBI = as.character(map$ChEBI)
  nodeAttrib.unmapp = nodeAttrib.unmapp %>%
    left_join(map[,c("ChEBI", "KEGG", "HMDB")], by = c("MappingID" = "HMDB")) %>%
    mutate(
      ChEBI = coalesce(as.character(ChEBI.x), as.character(ChEBI.y)),
      KEGG = coalesce(KEGG.x, KEGG.y)
    ) %>%
    dplyr::select(-ChEBI.x, -ChEBI.y, -KEGG.x, -KEGG.y)

  nodeAttrib$ChEBI = as.character(nodeAttrib$ChEBI)
  nodeAttrib = bind_rows(nodeAttrib[!x,], nodeAttrib.unmapp)

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
  m1 = match(network2$mappingID1, nodeAttrib$MappingID)
  m2 = match(network2$mappingID2, nodeAttrib$MappingID)

  diff_cols = grep("_diff", names(nodeAttrib), value = TRUE)

  for (col in diff_cols) {
    network2[[col]] = ifelse(
      nodeAttrib[[col]][m1] == "*" | nodeAttrib[[col]][m2] == "*",
      "*", ""
    )
  }
  edgeAttrib = network2
  edgeAttrib_solid = network2 %>% filter(EdgeLine == "Solid")
  print("Edge attribute file")

  return(list(
    nodeAttrib = nodeAttrib,
    edgeAttrib = edgeAttrib,
    edgeAttrib_solid = edgeAttrib_solid
  ))
}
