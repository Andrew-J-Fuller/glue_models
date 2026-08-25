# Amino Acid Acquirer 
AAA is a small script which takes a folder of CSV files exported from PRISM and returns a new CSV containing:
#### Column 1
  The protein of interest.
#### Column 2
  The amino acid sequence of the E3 ligase and it's HALO tag. 
#### Column 3
  The amino acid sequence of the NanoLuc tag and the protein of interest.
#### Column 4
  The resulting NANOBret signal ratio averaged over 8 replicates.

This formats the data for use in computational models, and automatically find the amino acid sequence for each assay tag and it's combination used (AA's acquired from UNIPROT)
