import pandas as pd

def find_repeats(path_to_db):
    df = pd.read_csv(path_to_db, sep="\t")
    boolean_series = df.duplicated(subset=["tax_id"])
    repeats = df["tax_id"][boolean_series]
    return repeats


def verify_db(path_to_db):
    try:
        repeats = find_repeats(path_to_db)
        if repeats.size > 0:
            raise ValueError
    except ValueError:
        repeated_id_str = ""
        for i in repeats:
            repeated_id_str += str(i)+", "
        print("There are one or more duplicate entries in this database:\n" + repeated_id_str[:-2])

verify_db("/Users/hubingbing/Desktop/emu/emu_database/taxonomy.tsv")

