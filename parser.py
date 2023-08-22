import os
import h5py
from uniprotein import UniProtein


class Parser:
    def __init__(self):
        self.proteins_file_name = "data/uniprotkb_AND_reviewed_true_2023_08_22.tsv"
        self.embeddings_file_name = "data/per-protein.h5"
        self.uni_proteins = {}

    # loads up the embedding from the h5 file for the given uniprotein entry
    # dependency on read_ecs to have been run first
    def read_embeddings(self, proportion=None):
        # UniProt is providing raw embeddings (per-protein and per-residue using the ProtT5 model)
        # for UniProtKB/Swiss-Prot and some reference proteomes of model organisms here:
        uniproteins_without_embeddings = 0
        with h5py.File(self.embeddings_file_name, "r") as file:
            file_size = os.path.getsize(file.name)
            print("opening embeddings file:" + file.name + ' with size:' + str(file_size))
            for entry, embedding in file.items():
                try:
                    protein = self.find_by_entry(entry)
                    protein.set_embedding(embedding)
                    self.uni_proteins[entry] = protein
                except KeyError:
                    print("no Uniprot for: " + entry)
                    uniproteins_without_embeddings += 1
                    pass
        print("number of uniproteins without embeddings: " + str(uniproteins_without_embeddings)
              + " out of total " + str(len(self.uni_proteins)))

    # creates UniProtein instances for the record in the tsv file
    def read_ecs(self):
        # UniProt is providing protein data here:
        # https://www.uniprot.org/uniprotkb?facets=reviewed%3Atrue&query=%2A
        # expecting a tsv file with the specified columns

        with open(self.proteins_file_name) as f:
            file_size = os.path.getsize(f.name)
            print("opening uniprot protein properties file:" + f.name + ' with size:' + str(file_size))
            header = next(f)
            for line in f:
                ws = line.split('\t')
                go_id = ws[9].strip().split(";")
                ecs_number = ws[8].strip().split(";")
                active_site = ws[7].strip().split(";")
                sequence = ws[6]
                length = ws[5]
                protein_name = ws[2]
                entry_name = ws[1]
                entry = ws[0]
                protein = UniProtein(entry, entry_name, protein_name,
                                     length, sequence, active_site, ecs_number, go_id)
                self.uni_proteins[entry] = protein

    def find_by_entry(self, entry):
        return self.uni_proteins[entry]

    def find_all(self):
        return self.uni_proteins
