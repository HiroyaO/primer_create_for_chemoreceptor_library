# Create primer for all receptors

# install primer3

pip install biopython primer3-py # opptional

# import libraries

import pandas as pd
from Bio.Seq import Seq
import primer3


# import the data

file_path = "input_genes.csv" # the table that include target genes
gene_column = "GeneSequence"  # the name of column of gene sequences

df = pd.read_csv(file_path) # you need to change if needed. Ex, elsx.

# Define the specific sequences to prepend to each primer

fwd_prefix = "NNNNNNNNNNNNNNNNNNNN"  # Specify the 20 bp prefix for forward primer
rev_prefix = "NNNNNNNNNNNNNNNNNNNN"  # Specify the 20 bp prefix for reverse primer


# Define a function to design primers for each gene sequence
def design_primers(sequence):

    # Get the length of the sequence
    sequence_length = len(sequence)
    
    # Primer design conditions
    primer_conditions = {
        'PRIMER_OPT_SIZE': 25,
        'PRIMER_MIN_SIZE': 17,
        'PRIMER_MAX_SIZE': 36,
        'PRIMER_OPT_TM': 63.0,
        'PRIMER_MIN_TM': 55.0,
        'PRIMER_MAX_TM': 66.0,
        'PRIMER_MIN_GC': 30,
        'PRIMER_MAX_GC': 70,
        'PRIMER_PRODUCT_SIZE_RANGE': [[sequence_length, sequence_length]],  # Amplify the entire sequence
    }
    
    # Design primers
    primer_results = primer3.bindings.designPrimers(
        {
            'SEQUENCE_TEMPLATE': str(sequence),
            'SEQUENCE_INCLUDED_REGION': [0, sequence_length]  # Target the entire sequence
        },
        primer_conditions
    )
    
    # Check if primer design was successful
    try:
        forward_primer = fwd_prefix + primer_results['PRIMER_LEFT_0_SEQUENCE']
        reverse_primer = rev_prefix + primer_results['PRIMER_RIGHT_0_SEQUENCE']
        forward_tm = primer_results['PRIMER_LEFT_0_TM']
        reverse_tm = primer_results['PRIMER_RIGHT_0_TM']
    except KeyError:
        # Return empty values if primer design failed
        forward_primer = reverse_primer = "Design Failed"
        forward_tm = reverse_tm = None
    
    return forward_primer, reverse_primer, forward_tm, reverse_tm

# Design primers and add the results as new columns
df[['Forward_Primer', 'Reverse_Primer', 'Forward_Tm', 'Reverse_Tm']] = df[gene_column].apply(
    lambda seq: pd.Series(design_primers(seq))
)

# Output the results as a new table
# Save the results as a new CSV file
output_file_path = "output_with_primers.csv"
df.to_csv(output_file_path, index=False)
