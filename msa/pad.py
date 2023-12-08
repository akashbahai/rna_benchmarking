from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys
import os

def pad_sequences(records):
    max_length = max(len(str(record.seq)) for record in records)
    padded_records = [record if len(str(record.seq)) == max_length else SeqRecord(
                            seq=Seq(str(record.seq).ljust(max_length, '-')),
                            id=record.id,
                            description=record.description)
                      for record in records]
    return padded_records

def write_fasta(output_file, records):
    with open(output_file, 'w') as file:
        SeqIO.write(records, file, 'fasta')

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python process_fasta_biopython.py <input_fasta_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    path, filename = os.path.split(input_file)
    filename_no_extension = filename.split('.')[0]
    output_file = os.path.join(path, filename_no_extension) + ".afa"

    records = list(SeqIO.parse(input_file, 'fasta'))
    padded_records = pad_sequences(records)
    write_fasta(output_file, padded_records)

    print(f"Processing complete. Padded sequences saved to {output_file}")
