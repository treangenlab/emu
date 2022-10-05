from os.path import exists
import pandas as pd


def find_repeats(path_to_tsv):
    df = pd.read_csv(path_to_tsv, sep="\t")
    boolean_series = df.duplicated(subset=["tax_id"])
    repeats = df["tax_id"][boolean_series]
    return set(repeats)


def verify_db(path_to_db):
    # Check whether both species_taxid.fasta and taxonomy.tsv are present
    if not (exists(path_to_db + "/" + "species_taxid.fasta") and exists(path_to_db + "/" + "taxonomy.tsv")):
        raise FileNotFoundError("Missing species_taxid.fasta or taxonomy.tsv")
    # Check for duplicates in taxonomy.tsv
    repeats = find_repeats(path_to_db + "/" + "taxonomy.tsv")
    if len(repeats) > 0:
        raise ValueError("There are one or more duplicate entries in this database:\n" + str(repeats))
    print("Database looks good!")

