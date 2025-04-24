import requests
import re
import csv
import time
from tqdm import tqdm

def fetch_kegg_reaction(reaction_id):
    url = f"http://rest.kegg.jp/get/{reaction_id}"
    try:
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        return response.text
    except requests.exceptions.RequestException as e:
        print(f"Failed to fetch reaction {reaction_id}: {e}")
        return None

def fetch_kegg_compound(compound_id):
    if not re.match(r'^C\d+$', compound_id):
        return 'Unknown'
    url = f"http://rest.kegg.jp/get/{compound_id}"
    try:
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        match = re.search(r'NAME\s+(.+?);', response.text)
        return match.group(1) if match else 'Unknown'
    except requests.exceptions.RequestException as e:
        print(f"Failed to fetch compound {compound_id}: {e}")
        return 'Unknown'

def parse_kegg_reaction(data):
    if not data:
        return []

    reaction_id = re.search(r'ENTRY\s+(\S+)', data).group(1)
    name_match = re.search(r'NAME\s+(.+)', data)
    equation_match = re.search(r'EQUATION\s+(.+)', data)
    enzyme_match = re.search(r'ENZYME\s+(.+)', data)

    if not equation_match:
        print(f"No equation found for reaction {reaction_id}")
        return []

    reaction_name = name_match.group(1) if name_match else 'Unknown'
    equation = equation_match.group(1)
    catalysts = enzyme_match.group(1).split() if enzyme_match else []

    def extract_compounds(equation_part):
        compounds = []
        for participant in equation_part.strip().split(' + '):
            participant = participant.strip()
            compound_id_match = re.search(r'(C\d+)', participant)
            if compound_id_match:
                compounds.append(compound_id_match.group(1))
        return compounds

    input_output = equation.split(' <=> ')
    inputs = extract_compounds(input_output[0])
    outputs = extract_compounds(input_output[1])

    parsed_data = []
    for participant_id in inputs:
        participant_name = fetch_kegg_compound(participant_id)
        parsed_data.append({
            "Reaction ID": reaction_id,
            "Reaction Name": reaction_name,
            "Participant ID": participant_id,
            "Participant Name": participant_name,
            "Type": "input"
        })
    for participant_id in outputs:
        participant_name = fetch_kegg_compound(participant_id)
        parsed_data.append({
            "Reaction ID": reaction_id,
            "Reaction Name": reaction_name,
            "Participant ID": participant_id,
            "Participant Name": participant_name,
            "Type": "output"
        })
    for catalyst in catalysts:
        parsed_data.append({
            "Reaction ID": reaction_id,
            "Reaction Name": reaction_name,
            "Participant ID": catalyst,
            "Participant Name": 'N/A',
            "Type": "catalyst"
        })

    return parsed_data

def save_to_csv(filename, data):
    with open(filename, 'w', newline='') as csvfile:
        fieldnames = ['Reaction ID', 'Reaction Name', 'Participant ID', 'Participant Name', 'Type']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for row in data:
            writer.writerow(row)

def process_reaction(reaction_id):
    data = fetch_kegg_reaction(reaction_id)
    if data is not None:
        return parse_kegg_reaction(data)
    return [
        {"Reaction ID": reaction_id, "Reaction Name": "Error", "Participant ID": "Error", "Participant Name": "Error",
         "Type": "Error"}]

def main():
    with open('inputs/reaction.txt') as file:
        reaction_ids = [line.split()[0] for line in file.readlines()]

    batch_size = len(reaction_ids) // 10
    for i in range(10):
        batch_ids = reaction_ids[i * batch_size: (i + 1) * batch_size]
        results = []
        try:
            for reaction_id in tqdm(batch_ids, desc=f"Processing batch {i + 1}/10"):
                result = process_reaction(reaction_id)
                results.extend(result)
                time.sleep(0.34)  # Add a delay to ensure no more than 3 requests per second
        except Exception as e:
            print(f"Error occurred during processing batch {i + 1}: {e}")
            break
        batch_filename = f'kegg_reactions_batch_{i + 1}.csv'
        save_to_csv(batch_filename, results)
        print(f"Batch {i + 1} complete. Results saved to '{batch_filename}'.")

    combined_results = []
    for i in range(10):
        batch_file = f'kegg_reactions_batch_{i + 1}.csv'
        with open(batch_file, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                combined_results.append(row)

    save_to_csv('outputs/kegg_reactions_combined.csv', combined_results)
    print("All batches processed and combined into 'kegg_reactions_combined.csv'.")

if __name__ == "__main__":
    main()
