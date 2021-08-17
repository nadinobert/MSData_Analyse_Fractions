# Analysis of MS data
# Collected fractions from SEC
# aim: check the abundance of OHR-related proteins from CBDB1 in the different fractions

import pandas as pd
import csv
import re
import os

# directory where result files are stored
directory_path = r'C:/Users/nadin/Documents/UFZ/2021 01 26 LB119a_SEC/20210208_MS_Data/'

# list of proteins to check for
ohr_proteins = ["rdhA", "omeB", "omeA", "rdhB", "hupX", "hupS", "hupL"]


def proteins_in_fractions(fraction, filename):
    # sucht in einem txt-file (filename) mit dem präfix "fraction" nach den ohr-proteinen aus cbdb1
    # gibt accession und abundnace des jeweiligen ohr proteins in einer liste (frac_result) aus
    # erster eintrage ist jeweils die fraction (präfix des filenames)
    proteins = csv.DictReader(open(directory_path + filename, "rt"), delimiter='\t')
    frac_result = {"fraction": fraction}
    frac_accession = []
    frac_protein = []
    frac_abundance = []
    pattern = re.compile(r'Abundance.*')  # different expression for abundance in every file! -> regex expression
    for row in proteins:  # row = dict
        protein_text = row.get("Description")
        for prot in ohr_proteins:
            if prot in protein_text:
                frac_accession.append(row.get("Accession"))
                frac_protein.append(prot)
                for (key) in row:
                    if pattern.fullmatch(key):
                        frac_abundance.append(row.get(key))
    frac_result["Accession"] = frac_accession
    frac_result["Protein"] = frac_protein
    frac_result["Abundance"] = frac_abundance
    return (frac_result)


def print_list_to_csv(result_list):
    with open(directory_path + "result.csv", 'w',
              newline='') as myfile:  # creates a new file in target directory with the name: "result" -> optional mit dem paramter "a" (creates a new file if possible. avoiding to überschreiben old file)
        wr = csv.writer(myfile)
        wr.writerows(result)
    # continue
# else:
# continue

def print_dict_to_csv(result_dict):
    with open(directory_path + "result.csv", 'w') as csv_file:
        csv_output = csv.writer(csv_file)
        csv_output.writerow(['fraction', 'Accession', 'Protein', 'Abundance'])

        for key in sorted(result_dict.keys()):
            csv_output.writerow([key] + result_dict[key])


def peptides_in_proteins(protein_list, filename):
    peptides = csv.DictReader(open(directory_path + filename, "rt"), delimiter='\t')
    data = list(peptides)
    result_list = []
    for elem in protein_list:
        peptide_list = []
        for row in data:
            if row.get("Master Protein Accessions") == elem:
                    peptide_list.append(row.get("Annotated Sequence"))
            else:
                continue
        result_list.append(peptide_list)
    return(result_list)

# collected results from all files
result = []

# iteration through all protein.txt and peptide.txt files in target directory
for file in os.listdir(directory_path):
    if file.endswith("_Proteins.txt"):
        frac = file.rsplit("_", 1)[0]  # setting the maxsplit paramater to 1, will return a list with two elements result -> z.B. "frac11_etcid"
        proteins_in_fraction = (proteins_in_fractions(frac, file))
        print(proteins_in_fraction)
        #result.append(proteins_in_fraction)

        #if len(
        #        proteins_in_fraction) > 1:  # wenn in der liste mehr als ein eintrag also mehr als z.B. frac11_etcid steht, suche die entsprechnden peptide der gefundene proteine
        #    peptide_file = frac + "_PeptideGroups.txt"
        #    peptides = (peptides_in_proteins(proteins_in_fraction, peptide_file))
        #    result.append(peptides)

#print(result)
print_list_to_csv(proteins_in_fraction)

