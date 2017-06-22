# Command to add column with source to GO files:

awk '{OFS="\t"} {print $0, "GO-0071305_cellular_response_to_vitamin_D.tsv"}' GO-0071305_cellular_response_to_vitamin_D.tsv > GO-0071305_cellular_response_to_vitamin_D_awk.tsv

# Then concatenated all:
cat *awk* > GO_VD_proteins.txt