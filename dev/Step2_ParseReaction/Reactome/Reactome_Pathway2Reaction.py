import requests
import pandas as pd
from tqdm import tqdm

# Load pathway list without column names
pathways_df = pd.read_csv('inputs/ReactomePathways.txt', sep='\t', header=None, names=['Pathway ID', 'Pathway Name', 'Species'])

# Filter pathways by species
species_of_interest = ['Homo sapiens', 'Mus musculus', 'Rattus norvegicus']
filtered_pathways_df = pathways_df[pathways_df['Species'].isin(species_of_interest)]

# Function to fetch reactions for a given pathway
def get_reactions_for_pathway(pathway_id):
    url = f"https://reactome.org/ContentService/data/pathway/{pathway_id}/containedEvents"
    try:
        response = requests.get(url)
        response.raise_for_status()
        return response.json()
    except requests.exceptions.HTTPError as e:
        print(f"HTTPError for pathway {pathway_id}: {e}")
        return []  # Return an empty list if the request fails
    except requests.exceptions.RequestException as e:
        print(f"RequestException for pathway {pathway_id}: {e}")
        return []  # Handle other request-related errors

# Initialize lists to hold the pathway-reaction relationships
data = {
    'Pathway ID': [],
    'Pathway Name': [],
    'Reaction ID': [],
    'Reaction stID': [],
    'Reaction stIDVersion': [],
    'Reaction Name': []
}

# Iterate through each pathway with a progress bar
for _, row in tqdm(filtered_pathways_df.iterrows(), total=filtered_pathways_df.shape[0], desc="Processing Pathways"):
    pathway_id = row['Pathway ID']
    pathway_name = row['Pathway Name']
    reactions = get_reactions_for_pathway(pathway_id)

    for reaction in reactions:
        try:
            reaction_id = reaction['dbId']
            reaction_stID = reaction['stId']
            reaction_stIDVersion = reaction['stIdVersion']
            reaction_name = reaction['displayName']

            data['Pathway ID'].append(pathway_id)
            data['Pathway Name'].append(pathway_name)
            data['Reaction ID'].append(reaction_id)
            data['Reaction stID'].append(reaction_stID)
            data['Reaction stIDVersion'].append(reaction_stIDVersion)
            data['Reaction Name'].append(reaction_name)
        except KeyError as e:
            print(f"KeyError for reaction in pathway {pathway_id}: {e}")

# Create a DataFrame from the data
result_df = pd.DataFrame(data)

# Save the DataFrame to a CSV file
result_df.to_csv('outputs/pathway_reaction_relationships.csv', index=False)