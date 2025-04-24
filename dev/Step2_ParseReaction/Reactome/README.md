# Reactome Parsing

This folder contains scripts and outputs for parsing Reactome pathway and reaction data for **human**, **mouse**, and **rat**.

------------------------------------------------------------------------

## Step 1: Pathway to Reaction

-   **Source:** Downloaded pathway list from [Reactome](https://reactome.org/download-data)
-   **Script:** `Reactome_Pathway2Reaction.py`
-   **Output:** `pathway_reaction_relationships.csv`

Maps pathways to their corresponding reactions.

------------------------------------------------------------------------

## Step 2: Parse Reactions

-   **Input:** `pathway_reaction_relationships.csv`
-   **Script:** `Reactome_ReactionDetails.py`
-   **Output:** `reaction_entities_details_v2.csv`

Parses inputs, outputs, and catalysts for each reaction.

------------------------------------------------------------------------

## Step 3: Parse Complexes

-   **Input:** `reaction_entities_details_v2.csv`
-   **Output:** `complex_participants_details.csv`

Recursively finds components of protein/compound complexes.

------------------------------------------------------------------------

## Step 4: Compile Results

-   **Script:** `Reactome_parsing.R`

Merges all outputs for downstream use.

------------------------------------------------------------------------

## Final outputs

-   `Reactome_Reactions_HSA.txt`
-   `Reactome_Reactions_MMU.txt`
-   `Reactome_Reactions_RNO.txt`
-   Filtered outputs for HSA, MMU, and RNO
