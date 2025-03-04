#This script corrects the terminal residues in the 3DRNA models

from Bio import PDB
import sys


def fix_residue_name(residue_name):
    standard_codes = {"ADE": "A", "CYT": "C", "GUA": "G", "URA": "U"}
    if residue_name in standard_codes:
        return standard_codes[residue_name]
    return residue_name

def remove_suffix(name):
    if name.endswith("3") or name.endswith("5"):
        print("Present in ",sys.argv[1])
        return name[1:-1]
    return name

def process_pdb(input_filename, output_filename):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("RNA", input_filename)

    for model in structure:
        for chain in model:
            for residue in chain:
                residue.resname = fix_residue_name(residue.resname)
                residue.resname = remove_suffix(residue.resname)

    io = PDB.PDBIO()
    io.set_structure(structure)
    io.save(output_filename)

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        print("Usage: python script.py input.pdb")
    else:
        input_filename = sys.argv[1]
        #output_filename = f"{input_filename.split('.')[0]}_corrected.pdb"
        output_filename = sys.argv[1]
        process_pdb(input_filename, output_filename)
        print("Residue names have been corrected and saved to", output_filename)