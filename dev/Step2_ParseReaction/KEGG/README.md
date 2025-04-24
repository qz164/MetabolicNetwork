# KEGG Parsing

This folder contains scripts to download and process KEGG reaction data, extract enzymes, and map genes.

------------------------------------------------------------------------

## Steps

### Step 1: Download Reactions

-   File: `reaction.txt`
-   Description: Reaction list from KEGG API

### Step 2: Parse Reactions

-   Script: `KEGG_reactions.py`
-   Extracts: Substrates, products, catalysts
-   Output: `kegg_reactions_combined.csv`
-   Note: Limit KEGG API calls to 3 per second

### Step 3: Get Enzymes and Genes

-   Source: [KEGG enzyme link](https://rest.kegg.jp/link/enzyme/hsa)
-   Script: `ECnumber2Gene_API.py`
-   Input: `enzyme_hsa_mmu.txt`
-   Output: gene symbols from KEGG API

### Step 4: ID Mapping & Filtering

-   Script: `kegg_reactions.R`
-   Maps compound names, filters for HSA/MMU/RNO

### Extra: UniProt Mapping (output file is used in step4)

-   Script: `ECnumber2Uniprot_Expasy.py`
-   Input: `enzyme.dat` from Expasy
-   Output: `enzyme_uniprot_mapping.csv`

------------------------------------------------------------------------

## Final outputs

-   `KEGG_Reactions_HSA.txt`
-   `KEGG_Reactions_MMU.txt`
-   `KEGG_Reactions_RNO.txt`
-   Filtered outputs for HSA, MMU, and RNO


