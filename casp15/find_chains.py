from Bio import PDB
import sys

def count_chains(pdb_file):
    # Create a PDB parser
    parser = PDB.PDBParser(QUIET=True)

    # Parse the PDB file
    structure = parser.get_structure("protein", pdb_file)

    # Get the list of models in the structure
    models = structure.get_models()

    # Initialize a set to store unique chain identifiers
    unique_chains = set()

    # Iterate over models and chains to count unique chains
    for model in models:
        chains = model.get_chains()
        for chain in chains:
            unique_chains.add(chain.get_id())

    return len(unique_chains)

if __name__ == "__main__":
    # Replace "your_file.pdb" with the path to your PDB file
    pdb_file_path = sys.argv[1]

    try:
        chain_count = count_chains(pdb_file_path)
        print(f"Number of chains in {pdb_file_path}: {chain_count}")
    except Exception as e:
        print(f"Error: {e}")