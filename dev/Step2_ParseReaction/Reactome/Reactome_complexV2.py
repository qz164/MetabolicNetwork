import requests
import pandas as pd
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed

# Load the list of complex IDs from the CSV file
reaction_entities_df = pd.read_csv('outputs/reaction_entities_details_v2.csv')

# Filter for the specified classes and determine the appropriate prefix based on the "Reaction ID"
specified_classes = ['DefinedSet', 'CandidateSet', 'Complex', 'ChemicalDrug', 'ProteinDrug', 'Polymer', 'RNADrug']
filtered_entries = reaction_entities_df[reaction_entities_df['Schema Class'].isin(specified_classes)][['Entity stID']]

# Get all specified stIDs
entity_ids = filtered_entries['Entity stID'].unique()

# Function to fetch detailed information for a given complex or reaction
def fetch_details(entity_id):
    url = f"https://reactome.org/ContentService/data/query/{entity_id}"
    retries = 3
    for _ in range(retries):
        try:
            response = requests.get(url, timeout=10)
            response.raise_for_status()
            return entity_id, response.json()
        except requests.exceptions.RequestException as e:
            err_msg = f"Error fetching details for entity ID {entity_id}: {e}"
            if _ == retries - 1:
                return entity_id, err_msg

# Initialize lists to hold the complex details
data = {
    'reaction_id': [],
    'reaction_name': [],
    'complex_id': [],
    'complex_name': [],
    'participant_id': [],
    'participant_name': [],
    'participant_type': [],
    'participant_stId': []
}

# Initialize list to hold errors
errors = []

# Number of concurrent threads
num_threads = 10

# Function to process participants and recursively fetch details if needed
def process_participants(participants, reaction_id, reaction_name, complex_id, complex_name):
    for participant in participants:
        if isinstance(participant, dict):
            schema_class = participant.get('schemaClass', '')
            if schema_class in specified_classes:
                # Recursively process components of the specified classes
                _, details = fetch_details(participant.get('stId'))
                if isinstance(details, str):  # Check if it is an error message
                    errors.append(details)
                    continue
                nested_participants = details.get('hasComponent', []) + \
                                      details.get('hasMember', []) + \
                                      details.get('hasCandidate', [])
                if nested_participants:
                    process_participants(nested_participants, reaction_id, reaction_name, complex_id, complex_name)
                else:
                    participant_id = participant.get('identifier', '') or participant.get('dbId', '')
                    participant_name = participant.get('displayName', '')
                    participant_type = participant.get('schemaClass', '')
                    participant_stId = participant.get('stId', '')

                    data['reaction_id'].append(reaction_id)
                    data['reaction_name'].append(reaction_name)
                    data['complex_id'].append(complex_id)
                    data['complex_name'].append(complex_name)
                    data['participant_id'].append(participant_id)
                    data['participant_name'].append(participant_name)
                    data['participant_type'].append(participant_type)
                    data['participant_stId'].append(participant_stId)
            else:
                participant_id = participant.get('identifier', '') or participant.get('dbId', '')
                participant_name = participant.get('displayName', '')
                participant_type = participant.get('schemaClass', '')
                participant_stId = participant.get('stId', '')

                data['reaction_id'].append(reaction_id)
                data['reaction_name'].append(reaction_name)
                data['complex_id'].append(complex_id)
                data['complex_name'].append(complex_name)
                data['participant_id'].append(participant_id)
                data['participant_name'].append(participant_name)
                data['participant_type'].append(participant_type)
                data['participant_stId'].append(participant_stId)

# Function to process a single entity ID
def process_entity(entity_id):
    # Fetch details for the entity
    _, entity_details = fetch_details(entity_id)

    if isinstance(entity_details, str):  # Check if it is an error message
        errors.append(entity_details)
        return

    entity_name = entity_details.get('displayName', '')

    # Process participants based on the schema class of the entity
    process_participants(entity_details.get('hasComponent', []) +
                         entity_details.get('hasMember', []) +
                         entity_details.get('hasCandidate', []), entity_id, entity_name, entity_id, entity_name)

# Process entities with a progress bar and concurrency
with ThreadPoolExecutor(max_workers=num_threads) as executor:
    futures = [executor.submit(process_entity, entity_id) for entity_id in entity_ids]
    for _ in tqdm(as_completed(futures), total=len(futures), desc="Processing Entities"):
        pass

# Create a DataFrame from the data
result_df = pd.DataFrame(data)

# Save the DataFrame to a CSV file with comma as the delimiter
result_df.to_csv('outputs/complex_participants_details.csv', index=False, sep=',')

# Save errors to a log file
with open('complex_details_errors.log', 'w') as error_file:
    for error in errors:
        error_file.write(f"{error}\n")

# Display the DataFrame
print(result_df)
