# R/merge_network.R

#' Merge Reactome, KEGG, and HMDB network and consolidate unique pairs
#'
#' @param reactome.edge Reactome edge attribute file
#' @param reactome.node Reactome node attribute file
#' @param kegg.edge Reactome edge attribute file
#' @param kegg.node Reactome node attribute file
#' @param hmdb.edge Reactome edge attribute file
#' @param hmdb.node Reactome node attribute file
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
merge_network <- function(reactome.edge,
                          reactome.node,
                          kegg.edge,
                          kegg.node,
                          hmdb.edge,
                          hmdb.node,
                          write_output = FALSE,
                          output_dir = "outputs") {

  # Use other internal datasets directly
  map = MetabolicNetwork:::map
  lipmap = MetabolicNetwork:::lipmap

  # Core processing
  result <- .merge_network_core(
    reactome.edge = reactome.edge,
    reactome.node = reactome.node,
    kegg.edge = kegg.edge,
    kegg.node = kegg.node,
    hmdb.edge = hmdb.edge,
    hmdb.node = hmdb.node,
    map = map,
    lipmap = lipmap
  )

  # Write outputs
  if (write_output) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

    utils::write.table(result$nodeAttrib,
                       file.path(output_dir, "Merged_nodes.txt"),
                       sep = "\t", quote = FALSE,
                       na = "", row.names = FALSE)

    utils::write.table(result$edgeAttrib,
                       file.path(output_dir, "Merged_edges.txt"),
                       sep = "\t", quote = FALSE,
                       na = "", row.names = FALSE)

    utils::write.table(result$edgeAttrib_solid,
                       file.path(output_dir, "Merged_edges_solid.txt"),
                       sep = "\t", quote = FALSE,
                       na = "", row.names = FALSE)
  }
}
.merge_network_core <- function(reactome.edge,
                                reactome.node,
                                kegg.edge,
                                kegg.node,
                                hmdb.edge,
                                hmdb.node,
                                map,
                                lipmap
    ){
  ##################################
  # 1. Metabolite ID harmonisation
  ##################################
  reactome.node = .unique_id(reactome.node, "Reactome")
  kegg.node = .unique_id(kegg.node, "KEGG")
  hmdb.node = .unique_id(hmdb.node, "HMDB")

  reactome.met = reactome.node[reactome.node$NodeType == "Metabolite",]
  kegg.met = kegg.node[kegg.node$NodeType == "Metabolite",]
  hmdb.met = hmdb.node[hmdb.node$NodeType == "Metabolite",]

  c = Reduce(intersect, list(names(reactome.met), names(kegg.met), names(hmdb.met))) # common columns
  all.met = rbind(reactome.met[,c], kegg.met[,c], hmdb.met[,c])
  all.met[is.na(all.met)] = ""
  tmp = unique(all.met$Name_map[all.met$detected == "*"])
  tmp = tmp[!is.na(tmp)]
  m = match(all.met$Name_map[all.met$detected == "*"],tmp)
  all.met$ID = ""
  all.met$ID[all.met$detected == "*"] = paste0("Exp_", m)

  x = map$ChEBI_ID_Charge
  x = unique(x[!is.na(x)])
  map = map[!map$ChEBI %in% x, ]

  map.tmp = map[map$KEGG %in% all.met$MappingID|
                  map$ChEBI %in% all.met$MappingID |
                  map$HMDB %in% all.met$MappingID |
                  map$ChEBI_ID_Charge %in% all.met$MappingID, ]

  map.tmp$HMDB[!map.tmp$HMDB %in% all.met$MappingID] = NA
  map.tmp$ChEBI[!map.tmp$ChEBI %in% all.met$MappingID] = NA
  map.tmp$KEGG[!map.tmp$KEGG %in% all.met$MappingID] = NA
  map.tmp$ChEBI_ID_Charge[!map.tmp$ChEBI_ID_Charge %in% all.met$MappingID] = NA

  x = map.tmp[,c("KEGG", "HMDB", "ChEBI", "ChEBI_ID_Charge")]
  x = !duplicated(x)
  map.tmp = map.tmp[x, c("Name", "KEGG", "HMDB", "ChEBI", "ChEBI_ID_Charge")]
  map.tmp = cbind("Map_No" = 1:nrow(map.tmp), map.tmp)
  map.tmp$ChEBI = as.character(map.tmp$ChEBI)

  map.tmp = map.tmp %>%
    mutate(across(c(HMDB, ChEBI, KEGG, ChEBI_ID_Charge), as.character))

  all.met = all.met %>%
    left_join(
      map.tmp %>%
        pivot_longer(
          cols = c(HMDB, ChEBI, KEGG, ChEBI_ID_Charge),
          values_drop_na = TRUE,
          values_to = "MappingID"
        ),
      by = "MappingID"
    )

  all.met = all.met%>%
    group_by(Map_No, ID) %>%
    mutate(across(
      where(~ all(grepl("^[-+]?[0-9]*\\.?[0-9]+$", na.omit(as.character(.))))),
      ~ as.numeric(as.character(.))
    )) %>%
    reframe(across(
      everything(),
      ~ {
        if (Map_No[1] == "" || is.na(Map_No[1])) return(as.character(.))
        v = if (is.character(.)) na_if(., "") else .
        paste(unique(na.omit(v)), collapse = ", ")
      }
    )) %>%
    ungroup() %>%
    mutate(across(everything(), ~ ifelse(is.na(.), "", .))) %>%
    mutate(Label_db = ifelse(grepl(",", NetworkID), Name, Label_db))

  all.met = .clean_mapping_data(all.met)
  all.met = all.met %>%
    group_by(Label_db, ID) %>%
    summarise(across(
      everything(),
      ~ {
        val = if (is.character(.)) na_if(., "") else .
        paste(unique(na.omit(val)), collapse = ", ")
      },
      .names = "{.col}"
    ), .groups = "drop") %>%
    mutate(across(everything(), ~ ifelse(is.na(.), "", .)))

  tmp = unlist(lapply(all.met$`HMDB lip`, function(x)strsplit(x, ", ")[[1]][1]))
  m = match(tmp, lipmap$accession)
  id = which(!is.na(m))
  all.met$`HMDB lip`[id] = lipmap$`HMDB lip`[m[!is.na(m)]]

  tmp = lapply(all.met$MappingID, function(x)unlist(strsplit(x, ", ")))
  tmp = sapply(tmp, function(x)x[grepl("HMDB_", x)])
  tmp = lapply(tmp, function(x) if (length(x) == 0) NA_character_ else x)
  tmp = unlist(tmp)
  id = which(!is.na(tmp))
  all.met$`HMDB lip`[id] = tmp[!is.na(tmp)]

  lip = all.met %>%
    filter(grepl("HMDB_", `HMDB lip`)) %>%
    group_by(`HMDB lip`) %>%
    filter(ID != "" | all(ID == "")) %>%
    summarise(across(
      everything(),
      ~ {
        val = if (is.character(.)) na_if(., "") else .
        paste(unique(na.omit(val)), collapse = ", ")
      },
      .names = "{.col}"
    ), .groups = "drop") %>%
    mutate(across(everything(), ~ ifelse(is.na(.), "", .)))

  all.met = all.met %>%
    filter(!grepl("HMDB_", `HMDB lip`)) %>%
    group_by(`HMDB lip`, ID, Name_map) %>%
    summarise(across(
      everything(),
      ~ {
        val = if (is.character(.)) na_if(., "") else .
        paste(unique(na.omit(val)), collapse = ", ")
      },
      .names = "{.col}"
    ), .groups = "drop") %>%
    mutate(across(everything(), ~ ifelse(is.na(.), "", .)))

  all.met = rbind(lip, all.met)
  all.met$Name_map = unlist(lapply(all.met$Name_map, function(x)strsplit(x, ", ")[[1]][1]))

  all.met = all.met %>%
    select(-ID, -Map_No, -Name, -name)

  ##################################
  # 2. Gene ID harmonisation
  ##################################
  reactome.prot = reactome.node[reactome.node$NodeType == "GeneProduct",]
  kegg.prot = kegg.node[kegg.node$NodeType == "GeneProduct",]
  hmdb.prot = hmdb.node[hmdb.node$NodeType == "GeneProduct",]

  c = Reduce(intersect, list(names(reactome.prot), names(kegg.prot), names(hmdb.prot))) # common columns
  all.prot = rbind(reactome.prot[,c], kegg.prot[,c], hmdb.prot[,c])
  all.prot[is.na(all.prot)] = ""

  all.prot = all.prot %>%
    group_by(Name_map) %>%
    reframe(
      across(where(is.character), ~ paste(unique(na.omit(na_if(., ""))), collapse = ", ")),
      across(where(Negate(is.character)), ~ paste(unique(na.omit(.)), collapse = ", "))
    ) %>%
    ungroup()
  all.prot = all.prot[!duplicated(all.prot),]

  ##################################
  # 3. Prepare node attribute file
  ##################################
  all.prot = as.data.frame(apply(all.prot, 2, as.character), stringsAsFactors = FALSE)
  all.met = as.data.frame(apply(all.met, 2, as.character), stringsAsFactors = FALSE)

  all.node = bind_rows(all.prot, all.met)
  all.node = cbind("NetworkID_All" = 1:nrow(all.node), all.node)

  nodeAttrib = all.node
  print("Node attribute file")

  ############################################################################
  ################### Reconstruct the network ################################
  ############################################################################
  mappingID = lapply(all.node$MappingID, function(x)unlist(strsplit(x, ", ")))

  # Precompute the mapping from MappingID to NetworkID_All
  mappingID_flat = unlist(mappingID)
  networkID_flat = rep(all.node$NetworkID_All, sapply(mappingID, length))

  # Apply the optimized replacement
  reactome.edge = .replace_ID(reactome.edge, mappingID_flat, networkID_flat)
  kegg.edge = .replace_ID(kegg.edge, mappingID_flat, networkID_flat)
  hmdb.edge = .replace_ID(hmdb.edge, mappingID_flat, networkID_flat)

  # update sorted_pairs
  reactome.edge = reactome.edge %>%
    mutate(sorted_pairs = paste0(pmin(mappingID1, mappingID2), "_", pmax(mappingID1, mappingID2))) %>%
    mutate(source = "Reactome")

  kegg.edge = kegg.edge %>%
    mutate(sorted_pairs = paste0(pmin(mappingID1, mappingID2), "_", pmax(mappingID1, mappingID2))) %>%
    mutate(source = "KEGG")

  hmdb.edge = hmdb.edge %>%
    mutate(sorted_pairs = paste0(pmin(mappingID1, mappingID2), "_", pmax(mappingID1, mappingID2))) %>%
    mutate(source = "HMDB")

  c = Reduce(intersect, list(names(reactome.edge), names(kegg.edge), names(hmdb.edge))) # common columns
  all.edge = rbind(reactome.edge[,c], kegg.edge[,c], hmdb.edge[,c])
  all.edge[is.na(all.edge)] = ""
  all.edge = all.edge[all.edge$mappingID1 != "" & all.edge$mappingID2 !="", ]


  # Define Relationships
  level_order = c("Enzyme", "Enzyme or transporter", "Transporter",
                  "Metabolic Transformation", "SameComplex", "Others", "Unknown")

  all.edge$mappingID1 = as.character(all.edge$mappingID1)
  all.edge$mappingID2 = as.character(all.edge$mappingID2)

  all.edge = all.edge %>%
    group_by(sorted_pairs) %>%
    summarise(across(everything(), ~ if (cur_column() == "EdgeColor") {
      highest = unique(.)  # Get unique values from the EdgeColor column
      # Match to the level_order and get the highest level (lowest numeric position)
      highest[which.min(match(highest, level_order, nomatch = length(level_order) + 1))]
    } else {
      paste(unique(na.omit(na_if(., ""))), collapse = ", ")
    })) %>%
    filter(if (all(grepl(", ", mappingID1))) row_number() == 1 else TRUE)

  firstID = function(xx){
    unlist(lapply(xx, function(x)strsplit(x, ", ")[[1]][1]))
  }

  all.edge$mappingID1 = firstID(all.edge$mappingID1)
  all.edge$mappingID2 = firstID(all.edge$mappingID2)
  all.edge$Label1 = all.node$Name_map[match(all.edge$mappingID1, all.node$NetworkID_All)]
  all.edge$Label2 = all.node$Name_map[match(all.edge$mappingID2, all.node$NetworkID_All)]
  all.edge[is.na(all.edge)] = ""

  m1 = match(all.edge$mappingID1, all.node$NetworkID_All)
  m2 = match(all.edge$mappingID2, all.node$NetworkID_All)

  diff_cols = grep("_diff", names(all.node), value = TRUE)
  for (col in diff_cols) {
    all.edge[[col]] = ifelse(
      all.node[[col]][m1] == "*" | all.node[[col]][m2] == "*",
      "*", ""
    )
  }
  edgeAttrib = all.edge
  edgeAttrib_solid = all.edge[all.edge$EdgeLine == "Solid", ]
  print("Edge attribute file")

  return(list(
    nodeAttrib = nodeAttrib,
    edgeAttrib = edgeAttrib,
    edgeAttrib_solid = edgeAttrib_solid
  ))
}
