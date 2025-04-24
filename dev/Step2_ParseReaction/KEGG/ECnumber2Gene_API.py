import requests
import pandas as pd
import time
from tqdm import tqdm

def get_gene_symbols_for_enzyme(ec_number):
    url = f"http://rest.kegg.jp/get/ec:{ec_number}"
    response = requests.get(url)
    if response.status_code != 200:
        raise Exception(f"Error querying KEGG API: {response.status_code}")
    return response.text

def extract_human_mouse_genes(data):
    gene_symbols = {"human": [], "mouse": [], "rat": []}
    lines = data.split("\n")
    in_genes_section = False
    for line in lines:
        if line.startswith("GENES"):
            in_genes_section = True
            if "HSA:" in line:
                human_genes = line.split("HSA:")[1].strip().split()
                gene_symbols["human"].extend(human_genes)
            if "MMU:" in line:
                mouse_genes = line.split("MMU:")[1].strip().split()
                gene_symbols["mouse"].extend(mouse_genes)
            if "RNO:" in line:
                rat_genes = line.split("RNO:")[1].strip().split()
                gene_symbols["rat"].extend(rat_genes)
            continue
        if in_genes_section:
            if line.startswith(" "):  # continuation line in the GENES section
                if "HSA:" in line:
                    human_genes = line.split("HSA:")[1].strip().split()
                    gene_symbols["human"].extend(human_genes)
                if "MMU:" in line:
                    mouse_genes = line.split("MMU:")[1].strip().split()
                    gene_symbols["mouse"].extend(mouse_genes)
                if "RNO:" in line:
                    rat_genes = line.split("RNO:")[1].strip().split()
                    gene_symbols["rat"].extend(rat_genes)
            elif not line.startswith(" "):  # end of the GENES section
                break
    return gene_symbols

def create_enzyme_gene_table(ec_numbers):
    records = []
    for idx, ec_number in enumerate(tqdm(ec_numbers, desc="Processing enzymes")):
        try:
            data = get_gene_symbols_for_enzyme(ec_number)
            genes = extract_human_mouse_genes(data)
            for species, symbols in genes.items():
                for symbol in symbols:
                    records.append({"enzyme": ec_number, "species": species, "gene_symbol": symbol})
            if (idx + 1) % 3 == 0:  # Limit to 3 requests per minute
                time.sleep(60)
        except Exception as e:
            print(f"Error processing EC number {ec_number}: {e}")

    if records:
        df = pd.DataFrame(records)
        df.to_csv("outputs/enzyme_gene_table.csv", index=False)
        print("Enzyme gene table saved to enzyme_gene_table.csv")
    else:
        print("No records to save.")

# Read the list of enzyme EC numbers from a text file
with open("inputs/enzyme_hsa_mmu_rno.txt", "r") as file:
    ec_numbers = [line.strip().split()[0] for line in file]

# Process the entire list of enzyme EC numbers
create_enzyme_gene_table(ec_numbers)
