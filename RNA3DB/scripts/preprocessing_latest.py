import os
from Bio.PDB import *
from Bio import SeqIO, PDB
import sys
from collections import defaultdict
import re



#this function returns a list of residues and the sequence of the pdb file
def get_residues(pdbfile):
    parser = PDBParser()
    structure = parser.get_structure('X', pdbfile)
    residues = []
    seq = ''
    mapping = {}
    for model in structure:
        for chain in model:
            for residue in chain:
                # Check if the residue is HETATM; if so, skip it
                # print(residue.get_id()[0].strip())
                if residue.get_id()[0].strip() == '':
                    res_id = str(residue.get_id()[1])+residue.get_id()[2].strip()
                    residues.append(res_id.strip())
                    seq += residue.get_resname()
                    mapping[res_id] = residue.get_resname()
    return seq, residues, mapping

#this function returns the list of residues that are missing in the pdb file
def get_missing_residues(pdbfile):
    h = parse_pdb_header(pdbfile)
    #print("Has missing residues:", h['missing_residues'])
    structure = PDBParser().get_structure('X', pdbfile)
    chains = [each.id for each in structure.get_chains()]

    #make a dictionary for each chain
    missing_per_chain = defaultdict(list)

    #go through each chain and get the missing residues
    for residue in h['missing_residues']:
        if residue["chain"] in chains:
            missing_per_chain[residue["chain"]].append(str(residue["ssseq"]).strip())
    
    for model in structure:
        for chain in model:
            for residue in chain:
                # Check if the residue is HETATM; if so, skip it
                # print(residue.get_id()[0].strip())
                hetres = ['H_GTP','H_TLN','H_LCA','H_LCG']
                if residue.get_id()[0].strip() in hetres:
                    missing_per_chain[chain.id].append(str(residue.get_id()[1]).strip())
    return chains, missing_per_chain

def format_fragments(residue_numbers):
    fragments = []
    current_fragment = []

    for num in residue_numbers:
        if not current_fragment or int(num) == int(current_fragment[-1]) + 1:
            current_fragment.append(num)
        else:
            fragments.append(current_fragment)
            current_fragment = [num]

    if current_fragment:
        fragments.append(current_fragment)

    formatted_fragments = []
    for fragment in fragments:
        if len(fragment) > 1:
            formatted_fragments.append(f"{fragment[0]}to{fragment[-1]}")
        else:
            formatted_fragments.append(fragment[0])

    return formatted_fragments

def custom_sort_key(s):
    parts = s.strip()
    
    # Use regular expressions to extract numeric and character parts
    numeric_part = re.search(r'-?\d+', parts).group(0)
    char_part = ''.join(filter(str.isalpha, parts))
    
    return (int(numeric_part), char_part)


def remap_pdb_residues(input_pdb_filename, mapping, real_seq):
    print(mapping)
    # Create a PDB parser
    pdb_parser = PDB.PDBParser(QUIET=True)

    # Load the input PDB structure
    structure = pdb_parser.get_structure("input", input_pdb_filename)
    new_structure = Structure.Structure("output")
    new_model = Model.Model(0)
    # Create a new structure
    # new_structure = Structure.Structure('output')
    chain_counter = 0
    new_chain = PDB.Chain.Chain('A')
    # Iterate through all residues in the structure
    for model in structure:
        for chain in model:
            if chain_counter > 0:
                break
            for i, residue in enumerate(list(chain)):
                
                # Check if the residue is HETATM; if so, skip it
                if residue.get_id()[0].strip() == '':

                    # Check if the residue number is within the new_residue_numbers list
                    if i < len(mapping.keys()):
                        # Update the residue ID with the new residue number
                        old_res_number = str(residue.get_id()[1])+residue.get_id()[2].strip()
                        new_res_number = mapping[old_res_number]
                        new_residue = residue.copy()

                        #print(i, residue.get_id()[1], mapping[residue.get_id()[1]])
                        #print(residue.get_resname(), real_seq[new_res_number-1])
                        if residue.get_resname() == real_seq[new_res_number-1]:
                            new_residue.id = (' ', new_res_number, ' ')
                            print(residue.id,new_residue.id)
                            new_chain.add(new_residue)
                            
                        else:
                            print("The residues are different!")
                            #print(residue.get_resname(), real_seq[new_res_number-1])
                            sys.exit()
                    else:
                        # If the list of new residue numbers is exhausted, break the loop
                        break
                    #print(new_residue)
                    
            chain_counter += 1
            new_model.add(new_chain)
            new_structure.add(new_model)

    # Save the modified structure to a new PDB file
    print(new_structure.child_list)
    io = PDB.PDBIO()
    io.set_structure(new_structure)
    io.save(input_pdb_filename.strip(".pdb") + "_remapped.pdb")


'''def remap_pdb_residues(input_pdb_filename, mapping, real_seq):
    # Create a PDB parser
    pdb_parser = PDB.PDBParser(QUIET=True)

    # Load the input PDB structure
    structure = pdb_parser.get_structure("input", input_pdb_filename)

    # Create a new structure
    new_structure = PDB.Structure.Structure("output")

    # Create a new chain
    new_chain = PDB.Chain.Chain('A')

    # Iterate through all residues in the structure
    for model in structure:
        for chain in model:
            for i, residue in enumerate(list(chain)):
                # Check if the residue is HETATM; if so, skip it
                if residue.get_id()[0].strip() == '':
                    # Check if the residue number is within the new_residue_numbers list
                    if i < len(mapping.keys()):
                        # Update the residue ID with the new residue number
                        old_res_number = str(residue.get_id()[1]) + residue.get_id()[2].strip()
                        new_res_number = mapping[old_res_number]

                        # Create a new residue with the updated ID
                        print(chain.id)
                        new_residue = PDB.Residue.Residue((' ', new_res_number, ' '), residue.get_resname(), 'A')

                        # Add atoms to the new residue (if applicable)
                        for atom in residue:
                            new_residue.add(atom.copy())

                        # Add the new residue to the new chain
                        new_chain.add(new_residue)

                    else:
                        # If the list of new residue numbers is exhausted, break the loop
                        break

    # Add the new chain to the new structure
    new_structure.add(new_chain)

    # Save the modified structure to a new PDB file
    io = PDB.PDBIO()
    io.set_structure(new_structure)
    output_filename = input_pdb_filename.strip(".pdb") + "_remapped.pdb"
    io.save(output_filename)'''

if __name__ == "__main__":
    pdbname = sys.argv[1] #pdb file
    fasta = sys.argv[2] #fasta file
    pdb_seq, res, pdb_mapping = get_residues(pdbname)
    chains, missing_per_chain = get_missing_residues(pdbname)

    if (len(chains)) > 1:
        print("More than one chain in pdb file. Please preprocess the file.")
        print("Only considering the first chain by default.")

    resindexes = sorted(list(set(sorted(res))), key=custom_sort_key)
    #print(pdb_seq, resindexes, pdb_mapping)

    missing_resindexes = sorted(list(set(missing_per_chain[chains[0]])), key=custom_sort_key)
    #print(missing_resindexes)

    for i in missing_resindexes:
        pdb_mapping[i] = '-'

    #print(pdb_mapping)

    #verify that the missing residues are not in the list of residues
    incorrect_missing = [el for el in missing_resindexes if el in resindexes]
    if len(incorrect_missing) > 0:
        print("The following missing residues are in the list of residues:", incorrect_missing)
        sys.exit()

    pdb_index = sorted(resindexes+missing_resindexes, key=custom_sort_key)
    #print(pdb_index)

    pdb_seq_complete = ''
    for i in pdb_index:
        pdb_seq_complete += pdb_mapping[i]
    #print(pdb_seq)
    #print(pdb_seq_complete)

    #here you should have your real sequence!
    real_seq = SeqIO.read(sys.argv[2], "fasta").seq

    fastamapping = defaultdict(list)
    fasta_index = []
    for i, j in enumerate(real_seq):
        fastamapping[i+1] = j
        fasta_index.append(i+1)

    #print(pdb_mapping, fastamapping)


    missed_seq = real_seq[:]
    #print(missed_seq)

    for i,el in enumerate(real_seq):
        if len(missing_resindexes) == 0:
            break
        #print(i, el)
        if pdb_index[i] in missing_resindexes:
            #print(i)
            missed_seq = missed_seq[:i]+'-'+missed_seq[i+1:]
    
    #print(missed_seq, "new")
    print(pdb_seq_complete)
    print(real_seq)

    if missed_seq  == pdb_seq_complete:
        print("The sequences are the same!")
    else:
        print("The sequences are different! Check",pdbname)
        

    #print(pdb_index, fasta_index)
    residue_mapping = {}
    if len(pdb_index) == len(fasta_index):
        for i in range(len(pdb_index)):
            residue_mapping[pdb_index[i]] = fasta_index[i]
    else:
        print("The sequences are different!")
        sys.exit()

    #print(residue_mapping)
    remap_pdb_residues(pdbname, residue_mapping, real_seq)

    pdb_res_seq = []
    for i in resindexes:
        pdb_res_seq.append(residue_mapping[i])

    native_fragments = format_fragments(pdb_res_seq)
    #print(native_fragments)
    model_fragments = []
    for fragment in native_fragments:
        line=fragment.split("to")
        if len(line) > 1:
            model_fragments.append(f"{str(line[0])}-{str(line[1])}")
    #print(model_fragments)
    
    fragment_file = pdbname.split("_native")[0] + "_fragments.txt"
    with open(fragment_file, "w") as file:
        file.write(f"\'A:")
        file.write('+'.join(model_fragments) + "\'")
    print("Renumbering and preprocessing done for:", pdbname)

