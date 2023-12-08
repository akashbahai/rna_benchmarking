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
                    residues.append(int(res_id.strip()))
                    seq += residue.get_resname()
                    mapping[int(res_id)] = residue.get_resname()
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
            missing_per_chain[residue["chain"]].append(int(str(residue["ssseq"]).strip()))
    
    for model in structure:
        for chain in model:
            for residue in chain:
                # Check if the residue is HETATM; if so, skip it
                # print(residue.get_id()[0].strip())
                res = ['A','U','G','C']
                if residue.get_id()[0].strip() in res:
                    missing_per_chain[chain.id].append(int(str(residue.get_id()[1]).strip()))
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
    parts = str(s).strip()
    
    # Use regular expressions to extract numeric and character parts
    numeric_part = re.search(r'-?\d+', parts).group(0)
    char_part = ''.join(filter(str.isalpha, parts))
    
    return (int(numeric_part), char_part)


def remap_pdb_residues(input_pdb_filename, mapping, real_seq):
    #print(mapping)
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
                    #print(list(mapping)[0])
                    if i < len(mapping.keys()) + int(list(mapping)[0]):
                        # Update the residue ID with the new residue number
                        old_res_number = str(residue.get_id()[1])+residue.get_id()[2].strip()
                        new_res_number = mapping[old_res_number]
                        new_residue = residue.copy()

                        #print(i, residue.get_id()[1], mapping[residue.get_id()[1]])
                        #print(residue.get_resname(), real_seq[new_res_number-1])
                        new_residue.id = (' ', new_res_number, ' ')
                        #print(residue.id,new_residue.id)
                        new_chain.add(new_residue)
                        
                    else:
                        # If the list of new residue numbers is exhausted, break the loop
                        break
                    #print(new_residue)
                    
            chain_counter += 1
            new_model.add(new_chain)
            new_structure.add(new_model)

    # Save the modified structure to a new PDB file
    #print(new_structure.child_list)
    io = PDB.PDBIO()
    io.set_structure(new_structure)
    #print(input_pdb_filename.strip(".pdb"))
    io.save(input_pdb_filename.split('_')[0] + "_native_remapped.pdb")



def find_starting_index(string1, string2):
    # Replace underscores in string2 with a regex pattern for any character
    pattern = re.sub(r'_', '.', string2)
    
    # Use re.search to find the starting index of the pattern in string1
    match = re.search(pattern, string1)
    
    # Return the starting index or -1 if not found
    return match.start() if match else -1

def largest_common_substring(string1, string2, numbering1, numbering2):
    len1 = len(string1)
    len2 = len(string2)

    # Initialize a 2D array to store lengths of common suffixes
    dp = [[0] * (len2 + 1) for _ in range(len1 + 1)]

    # Variables to store length of the longest common substring and its ending position
    max_len = 0
    end_pos1 = 0
    end_pos2 = 0

    for i in range(1, len1 + 1):
        for j in range(1, len2 + 1):
            if string1[i - 1] == string2[j - 1] or string2[j - 1] == '_':
                dp[i][j] = dp[i - 1][j - 1] + 1
                if dp[i][j] > max_len:
                    max_len = dp[i][j]
                    end_pos1 = i
                    end_pos2 = j
            else:
                dp[i][j] = 0

    # Extract the largest common substring from string2
    common_substring = string2[end_pos2 - max_len : end_pos2]

    # Extract the numbering of the common substring
    numbering_common_substring = numbering2[end_pos2 - max_len : end_pos2]

    # Extract the numbering of the starting position in both strings
    start_numbering1 = numbering1[end_pos1 - max_len]
    start_numbering2 = numbering2[end_pos2 - max_len]

    return (
        common_substring,
        numbering_common_substring,
        start_numbering1,
        start_numbering2
    )

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
        pdb_mapping[i] = '_'

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


    common=''
    index_common=[]
    residue_mapping={}
    if missed_seq  == pdb_seq_complete:
        print("The sequences are the same!")
        common = missed_seq
        index_common = pdb_index[0:len(pdb_index)]
        for i in range(len(pdb_index)):
                residue_mapping[str(pdb_index[i])] = fasta_index[i]
    else:
        print("The sequences are different! Look at residue numbering as well. Check",pdbname)
        testseq=''
        #print(pdb_mapping)
        for i in range(pdb_index[0],pdb_index[-1]+1):
            try:
                testseq += pdb_mapping[i]
            except:
                testseq += '_'
                pdb_mapping[i] = '_'
                missing_resindexes.append(int(i))
        print(testseq)
        print(real_seq)
        #print(fasta_index)
        pdb_index = sorted(list(set(pdb_index) | set (missing_resindexes)))

        common, index_common, startingindex_1, startingindex_2 = largest_common_substring(real_seq, testseq, fasta_index, pdb_index)
        #print(real_seq)
        #print(testseq)
        print(common)
        #print(startingindex_1, startingindex_2, common)
        startinggaps = startingindex_1 - 1
        endinggaps = len(real_seq) - len(common) - startingindex_1 + 1
        aligned_testseq = '_' * (startinggaps) + common + '_' * (endinggaps)
        print(real_seq)
        print(aligned_testseq)
        print(startingindex_1, startingindex_2)

        

        #update the missing_resindexes
        # Extend on the left
        print(missing_resindexes)
        print(index_common)
        if len(missing_resindexes)>0:
            missing_resindexes[:0] = range(pdb_index[0], pdb_index[0] + startinggaps - (len(common) - len(common.lstrip('_'))))
        else:
            missing_resindexes = list(range(pdb_index[0], pdb_index[0] + startinggaps - (len(common) - len(common.lstrip('_')))))

        # Extend on the right
        missing_resindexes.extend(range(index_common[-1]+1, pdb_index[-1]+1))
        print(missing_resindexes)
        missing_resindexes = sorted(list(set(missing_resindexes)))

        pdb_index_pseudo = sorted(list(set(pdb_index) | set (missing_resindexes)))
        print(pdb_index_pseudo)
        shift_residues_by = startingindex_1 - startingindex_2
        print(shift_residues_by, startingindex_2)
        #print(pdb_index, fasta_index)
        residue_mapping = {}
        for i in range(len(pdb_index_pseudo)):
                residue_mapping[str(pdb_index_pseudo[i])] = pdb_index_pseudo[i] + shift_residues_by
    

    print(residue_mapping)
    remap_pdb_residues(pdbname, residue_mapping, real_seq)

    common_residues = []
    print(index_common, len(index_common), len(common), common)
    for i in range(0,len(common)):
        if common[i]!= '_':
            #print(common[i])
            common_residues.append(index_common[i])
    print(common_residues)

    pdb_res_seq = []
    for i in common_residues:
        pdb_res_seq.append(int(residue_mapping[str(i)]))

    native_fragments = format_fragments(pdb_res_seq)
    #print(native_fragments)
    model_fragments = []
    for fragment in native_fragments:
        print(fragment)
        try:
            line=fragment.split("to")
            if len(line) > 1:
                model_fragments.append(f"{str(line[0])}-{str(line[1])}")
        except:
            model_fragments.append(str(fragment))
        #print(model_fragments)
    
    fragment_file = pdbname.split("_native")[0] + "_fragments.txt"
    with open(fragment_file, "w") as file:
        file.write(f"\'A:")
        file.write('+'.join(model_fragments) + "\'")
    print("Renumbering and preprocessing done for:", pdbname)

