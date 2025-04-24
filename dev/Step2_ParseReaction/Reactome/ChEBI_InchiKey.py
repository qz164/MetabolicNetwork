# /Users/zhaoqing/Library/CloudStorage/Dropbox/Macbook/DataResource/Reactome/pythonProject/ChEBI_InchiKey.py

import os


def parse_sdf(sdf_file):
    chebi_data = []
    with open(sdf_file, 'r') as f:
        lines = f.readlines()

    mol_block = []
    in_mol_block = False
    for line in lines:
        if line.startswith('$$$$'):
            chebi_id = None
            chebi_name = None
            inchi_key = None
            kegg_links = []
            other_data = {}

            for mol_line in mol_block:
                if mol_line.startswith('> <ChEBI ID>'):
                    chebi_id_index = mol_block.index(mol_line) + 1
                    chebi_id = mol_block[chebi_id_index].strip()
                elif mol_line.startswith('> <ChEBI Name>'):
                    chebi_name_index = mol_block.index(mol_line) + 1
                    chebi_name = mol_block[chebi_name_index].strip()
                elif mol_line.startswith('> <InChIKey>'):
                    inchi_key_index = mol_block.index(mol_line) + 1
                    inchi_key = mol_block[inchi_key_index].strip()
                elif mol_line.startswith('> <KEGG COMPOUND Database Links>'):
                    kegg_index = mol_block.index(mol_line) + 1
                    kegg_links.append(mol_block[kegg_index].strip())
                elif mol_line.startswith('> <'):
                    tag = mol_line.strip()[3:-1]
                    data_index = mol_block.index(mol_line) + 1
                    other_data[tag] = mol_block[data_index].strip()

            if chebi_id:
                chebi_data.append({
                    "ChEBI_ID": chebi_id,
                    "ChEBI_Name": chebi_name,
                    "InChIKey": inchi_key,
                    "KEGG_Links": ", ".join(kegg_links),  # Join list items with a comma
                    **other_data
                })
            mol_block = []
            in_mol_block = False
        else:
            mol_block.append(line)
            if not in_mol_block and line.strip() == '':
                in_mol_block = True

    return chebi_data


def write_to_tsv(data, output_file):
    headers = set()
    for entry in data:
        headers.update(entry.keys())

    headers = sorted(headers)

    with open(output_file, 'w') as f:
        f.write("\t".join(headers) + "\n")
        for entry in data:
            f.write("\t".join(str(entry.get(header, "")) for header in headers) + "\n")


def main():
    sdf_file = 'ChEBI_complete_3star.sdf'
    output_file = 'ChEBI_data.tsv'

    if not os.path.exists(sdf_file):
        print(f"File not found: {sdf_file}")
        return

    chebi_data = parse_sdf(sdf_file)
    write_to_tsv(chebi_data, output_file)
    print(f"Output written to {output_file}")


if __name__ == '__main__':
    main()
