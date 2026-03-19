# RABV Variant Caller: Protein-Level Hierarchical Typing

A Python pipeline for codon-aware variant calling and hierarchical phylogenetic typing from viral sequencing data.

## Key Features
- **Codon-aware variant calling**: Identifies amino acid changes at protein level
- **Hierarchical typing**: Groups samples based on maximal shared mutation patterns
- **Multi-level classification**: Generates lineage assignments (A → A.1 → A.1.2)
- **Reference flexible**: Works with any virus given a GenBank annotation file
- **Parallel processing**: Multi-threaded for large datasets
- **Comprehensive output**: Excel file with mutation summaries and hierarchical types

## Installation
```bash
git clone https://github.com/yourusername/rabv-variant-caller.git
cd rabv-variant-caller
pip install -r requirements.txt


## Quick Start
python rabv_variant_caller.py -g reference.gb -r reference.fasta -s samples.fasta -o output

Input Requirements

GenBank file (.gb): CDS annotations with coordinates

    Reference FASTA: Reference genome sequence

    Sample FASTA: Sequences to analyze

## Sheet	Description
Type_Definitions	Hierarchical types with level-specific mutations
Mutation_Levels	Every mutation classified by level
Sample_Details	Per-sample mutation presence
Group_Hierarchy	Hierarchical group structure
Mutation_Frequencies	All mutations with frequencies
Protein_Summary	Mutation statistics by protein
Parameters	Analysis parameters used


## Hierarchical Typing Logic

 * Level 1 (Major Groups): Maximal mutations shared by ALL samples in group
   * Largest group = 'A', next largest = 'B', etc.
   * Mutations appear ONLY in defining_mutations column
 * Level 2 (Subgroups): NEW mutations within each major group
   * Simpler groups = '.1', more complex = '.2', '.3'
   * Mutations NEVER appear in higher levels
 * Level 3+: Continue recursively until all mutations used


```bash
Example:
A (defining: G:V13A, G:F14S)
├── A.1 (defining only)
└── A.2 (additional: N:N1838D, L:P2121S)
    ├── A.2.1 (defining only)
    └── A.2.2 (additional mutations)
```

|  **Parameter**        |	          **Description**     |     **Default**     |
|-----------------------|-------------------------------|---------------------|
|--max-levels	          |  Maximum hierarchical levels	|          5          |
|--threshold	          |  Conserved mutation frequency |        	0.9         |
|--singleton-threshold  |	Min frequency for typing      |	       0.05         |
|--min-group-size	      | Min samples per group         |	        2           |
|--include-conserved	  |  Include conserved mutations  |	       False        |
|-t, --threads	        | Number of threads	            |      CPU count      |

## Type Definition Columns
|         **Column**	|         **Description**                    |
|-----------------------|--------------------------------------------|
|type	                | Hierarchical designation (A, A.1, B.2.1)   |
|Level1, 2, 3	        | Individual level components                |
|samples	            |    Samples in this type                    |
|count	                |    Number of samples                       |
|defining_mutations	    |     Level 1 mutations                      |
|level2/3_mutations	    |    Subgroup mutations                      |
|total_mutations        |	Combined mutation count                  |

## Dependencies
|Package  |	Version  |
|---------|----------|
|Python	  | ≥ 3.8    |
|biopython|	≥ 1.79   |
|pandas	  | ≥ 1.3.0  |
|openpyxl	| ≥ 3.0.9  |
|tqdm	    | ≥ 4.62.0 |
