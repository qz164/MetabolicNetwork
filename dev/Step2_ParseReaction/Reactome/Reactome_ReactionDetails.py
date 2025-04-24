import requests
import pandas as pd
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed

# Load the list of reaction IDs from the CSV file
pathway_reaction_df = pd.read_csv('outputs/pathway_reaction_relationships.csv')
reaction_ids = pathway_reaction_df['Reaction stID'].unique()

# Function to fetch detailed information for a given reaction with retries
def get_reaction_details(reaction_id):
    url = f"https://reactome.org/ContentService/data/query/{reaction_id}"
    retries = 3
    for _ in range(retries):
        try:
            response = requests.get(url, timeout=10)
            response.raise_for_status()
            return reaction_id, response.json()
        except requests.exceptions.RequestException as e:
            err_msg = f"Error fetching details for reaction ID {reaction_id}: {e}"
            if _ == retries - 1:
                return reaction_id, err_msg

# Function to fetch details for a given entity with retries
def get_entity_details(entity_id):
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

# Initialize lists to hold the reaction details
data = {
    'Reaction ID': [],
    'Reaction Name': [],
    'Class Name': [],
    'Entity ID': [],
    'Entity Name': [],
    'Schema Class': [],
    'Type': [],
    'Entity stID': []
}

# Initialize list to hold errors
errors = []

# Number of concurrent threads
num_threads = 10

# Process reactions with a progress bar
with ThreadPoolExecutor(max_workers=num_threads) as executor:
    future_to_reaction_id = {executor.submit(get_reaction_details, reaction_id): reaction_id for reaction_id in reaction_ids}
    for future in tqdm(as_completed(future_to_reaction_id), total=len(reaction_ids), desc="Processing Reactions"):
        reaction_id, reaction_details = future.result()

        if isinstance(reaction_details, str):  # Check if it is an error message
            errors.append(reaction_details)
            continue

        reaction_name = reaction_details.get('displayName', '')
        class_name = reaction_details.get('className', '')

        def process_entity(entity, entity_type):
            if isinstance(entity, dict):
                entity_id = entity.get('dbId', '')
                entity_stID = entity.get('stId', '')
                entity_name = entity.get('displayName', '')
                schema_class = entity.get('schemaClass', '')

                data['Reaction ID'].append(reaction_id)
                data['Reaction Name'].append(reaction_name)
                data['Class Name'].append(class_name)
                data['Entity ID'].append(entity_id)
                data['Entity Name'].append(entity_name)
                data['Schema Class'].append(schema_class)
                data['Type'].append(entity_type)
                data['Entity stID'].append(entity_stID)
            else:
                errors.append(f"Unexpected entity format: {entity}")

        # Process inputs
        for input_entity in reaction_details.get('input', []):
            if isinstance(input_entity, dict):
                process_entity(input_entity, 'input')

        # Process outputs
        for output_entity in reaction_details.get('output', []):
            if isinstance(output_entity, dict):
                process_entity(output_entity, 'output')

        # Process catalysts
        for catalyst in reaction_details.get('catalystActivity', []):
            if isinstance(catalyst, dict):
                catalyst_id = catalyst.get('dbId', '')
                _, catalyst_details = get_entity_details(catalyst_id)
                if isinstance(catalyst_details, str):  # Check if it is an error message
                    errors.append(catalyst_details)
                    continue

                physical_entity = catalyst_details.get('physicalEntity')
                if physical_entity:
                    if isinstance(physical_entity, dict):
                        process_entity(physical_entity, 'catalyst')
                        for sub_entity in physical_entity.get('hasComponent', []) + \
                                         physical_entity.get('hasMember', []) + \
                                         physical_entity.get('hasCandidate', []):
                            process_entity(sub_entity, 'catalyst')
                    elif isinstance(physical_entity, list):
                        for entity in physical_entity:
                            process_entity(entity, 'catalyst')
                            for sub_entity in entity.get('hasComponent', []) + \
                                             entity.get('hasMember', []) + \
                                             entity.get('hasCandidate', []):
                                process_entity(sub_entity, 'catalyst')

        # Process regulators
        for regulator in reaction_details.get('regulator', []):
            if isinstance(regulator, dict):
                process_entity(regulator, 'regulator')

        # Process active units if they exist
        for active_unit in reaction_details.get('activeUnit', []):
            if isinstance(active_unit, dict):
                process_entity(active_unit, 'activeUnit')

        # Process physical entities if they exist
        for physical_entity in reaction_details.get('physicalEntity', []):
            if isinstance(physical_entity, dict):
                process_entity(physical_entity, 'physicalEntity')

        # Process drugs if they exist
        for drug in reaction_details.get('drug', []):
            if isinstance(drug, dict):
                process_entity(drug, 'drug')

# Create a DataFrame from the data
result_df = pd.DataFrame(data)

# Save the DataFrame to a CSV file
result_df.to_csv('outputs/reaction_entities_details_v2.csv', index=False)

# Save errors to a log file
with open('reaction_details_errors_v2.log', 'w') as error_file:
    for error in errors:
        error_file.write(f"{error}\n")

# Display the DataFrame
print(result_df)
