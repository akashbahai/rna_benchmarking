import sys

def replace_t_with_u(seq):
    # Replace 'T' with 'U' in a DNA sequence
    return seq.replace('T', 'U')

if len(sys.argv) != 3:
    print("Usage: python test.py input_file.a3m output_file.a3m")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

# Read the a3m formatted file
with open(input_file, 'r') as file:
    lines = file.readlines()

# Modify DNA sequences by replacing 'T' with 'U'
modified_lines = []
for line in lines:
    if line.startswith(">"):
        modified_lines.append(line)
    else:
        modified_lines.append(replace_t_with_u(line))

# Write the modified alignment to a new a3m file
with open(output_file, 'w') as file:
    file.writelines(modified_lines)

print("Modified alignment saved to", output_file)
