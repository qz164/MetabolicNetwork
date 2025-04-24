# R/utils.R

#' Recursively find ontology offspring
#'
#' @keywords internal
#' @noRd
.find_offspring = function(chebi_id, df, visited = NULL) {
  if (is.null(visited)) {
    visited = list()
  }
  # Add the current chebi_id to the visited list
  visited = c(visited, chebi_id)

  children = df %>%
    filter(ID2 == chebi_id) %>%
    pull(ID1) # ID1 : child of ID2
  # Recursively find offspring for each child
  for (child in children) {
    if (!(child %in% visited)) {
      visited = .find_offspring(child, df, visited)
    }
  }
  return(visited)
}

#' Concatenate values and keep unique
#'
#' @keywords internal
#' @noRd
.concat_unique = function(x) {
  x = unique(x[!is.na(x) & x != ""])
  paste(x, collapse = ", ")
}

#' Customize network size
#'
#' @keywords internal
#' @noRd
.clean_metabolic_network_reactome = function(network_data,
                                    type = c("MetMN_Diff", "MetMN_Expand", "MN_Diff")) {
  type = match.arg(type)
  network3 = switch(type,
                    "MetMN_Diff" = {
                      network_data %>%
                        dplyr::filter(SigNode1 == "*" & SigNode2 == "*") %>%
                        dplyr::ungroup() %>%
                        dplyr::filter(!(grepl("[A-Za-z]", mappingID1) & grepl("[A-Za-z]", mappingID2)))
                    },
                    "MetMN_Expand" = {
                      # Get seed nodes from MetMN_Diff criteria
                      seed_network = network_data %>%
                        dplyr::filter(SigNode1 == "*" & SigNode2 == "*") %>%
                        dplyr::ungroup() %>%
                        dplyr::filter(!(grepl("[A-Za-z]", mappingID1) & grepl("[A-Za-z]", mappingID2)))
                      # Get all unique significant IDs
                      sig_ids = unique(c(seed_network$mappingID1, seed_network$mappingID2))
                      # Expand network to include neighbors
                      network_data %>%
                        dplyr::filter(mappingID1 %in% sig_ids | mappingID2 %in% sig_ids)
                    },
                    "MN_Diff" = {
                      network_data %>%
                        dplyr::filter(SigNode1 == "*" & SigNode2 == "*") %>%
                        dplyr::ungroup()
                    }
  )
  # Remove duplicate pairs and same complex relationships
  network3 %>%
    dplyr::group_by(`Reaction ID`, sorted_pairs) %>%
    dplyr::filter(!(EdgeColor == "" & any(EdgeColor == "SameComplex"))) %>%
    dplyr::slice(if (all(EdgeColor == "SameComplex")) 1 else dplyr::row_number()) %>%
    dplyr::ungroup() %>%
    return()
}

#' Customize network size
#'
#' @keywords internal
#' @noRd
.clean_metabolic_network_kegg = function(network_data,
                                   type = c("MetMN_Diff", "MetMN_Expand", "MN_Diff")) {
  type = match.arg(type)
  network3 = switch(type,
                    "MetMN_Diff" = {
                      network_data %>%
                        dplyr::filter(SigNode1 == "*" & SigNode2 == "*") %>%
                        dplyr::ungroup() %>%
                        dplyr::filter(!grepl("^C\\d{5}$", mappingID1) & !grepl("^C\\d{5}$", mappingID2))
                    },
                    "MetMN_Expand" = {
                      # Get seed nodes from MetMN_Diff criteria
                      seed_network = network_data %>%
                        dplyr::filter(SigNode1 == "*" & SigNode2 == "*") %>%
                        dplyr::ungroup() %>%
                        dplyr::filter((grepl("^C\\d{5}$", mappingID1) | grepl("^C\\d{5}$", mappingID2)))
                      # Get all unique significant IDs
                      sig_ids = unique(c(seed_network$mappingID1, seed_network$mappingID2))
                      # Expand network to include neighbors
                      network_data %>%
                        dplyr::filter(mappingID1 %in% sig_ids | mappingID2 %in% sig_ids)
                    },
                    "MN_Diff" = {
                      network_data %>%
                        dplyr::filter(SigNode1 == "*" & SigNode2 == "*") %>%
                        dplyr::ungroup()
                    }
  )
  # Remove duplicate pairs and same complex relationships
  network3 %>%
    dplyr::group_by(`Reaction ID`, sorted_pairs) %>%
    dplyr::ungroup() %>%
    return()
}

#' Clean mapping data for network merging
#'
#' @keywords internal
#' @noRd
.clean_mapping_data = function(data) {
  data = data %>%
    mutate(SplitMappingID = strsplit(MappingID, ",\\s*"))  # split by comma and possible spaces

  is_subset = function(list1, list2) {
    all(list1 %in% list2)}

  to_remove = logical(nrow(data))
  for (i in 1:nrow(data)) {
    for (j in 1:nrow(data)) {
      if (i != j && is_subset(data$SplitMappingID[[i]], data$SplitMappingID[[j]])) {
        to_remove[i] = TRUE
        break  # if it's a subset, no need to check further for this row
      }
    }
  }
  data_cleaned = data[!to_remove, ] %>% select(-SplitMappingID)
  return(data_cleaned)
}


#' Create a lookup function
#'
#' @keywords internal
#' @noRd
.replace_id_vectorized = function(id, mappingID_flat, networkID_flat) {
  matched_index = match(TRUE, mappingID_flat %in% id)
  if (!is.na(matched_index)) {
    return(networkID_flat[matched_index])
  } else {
    return(NA)  # Return NA if no match is found
  }
}

#'Vectorized function to replace IDs in the network
#'
#' @keywords internal
#' @noRd
.replace_ID = function(network, mappingID_flat, networkID_flat) {
  network$mappingID1 = sapply(network$mappingID1, .replace_id_vectorized, mappingID_flat, networkID_flat)
  network$mappingID2 = sapply(network$mappingID2, .replace_id_vectorized, mappingID_flat, networkID_flat)
  # Need to solve the problem of NA in ID mapping --> for now we delete those "pairs"
  network = network[!is.na(network$mappingID1)&!is.na(network$mappingID2), ]
  return(network)
}



#' Assign unique IDs for network merging
#'
#' @keywords internal
#' @noRd
.unique_id = function(node.data, db){
  if (db == "Reactome"){
    i = "R"
  } else if (db == "KEGG"){
    i = "K"
  } else if (db == "HMDB"){
    i = "H"
  }
  x = paste0(rep(i, nrow(node.data)), 1:nrow(node.data))
  node.data = cbind("NetworkID" = x, node.data)
  return(node.data)
}
