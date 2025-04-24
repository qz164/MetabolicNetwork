import re
import csv

def parse_enzyme_dat(file_path):
    ec_number = None
    entries = []
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('ID'):
                ec_number = line.split()[1].strip()
            elif line.startswith('DR') and ec_number:
                dr_entries = line.strip().split('DR   ')[1].split(';')
                for entry in dr_entries:
                    entry = entry.strip()
                    if entry:
                        uniprot_id, organism = entry.split(', ')
                        entries.append({
                            'EC Number': ec_number,
                            'UniProt ID': uniprot_id.strip(),
                            'Organism': organism.split('_')[1]
                        })
            elif line.startswith('//'):
                ec_number = None
    return entries

def save_to_csv(filename, data):
    with open(filename, 'w', newline='') as csvfile:
        fieldnames = ['EC Number', 'UniProt ID', 'Organism']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for row in data:
            writer.writerow(row)

def main():
    input_file = 'enzyme.dat'
    output_file = 'enzyme_uniprot_mapping.csv'
    data = parse_enzyme_dat(input_file)
    save_to_csv(output_file, data)

if __name__ == "__main__":
    main()
