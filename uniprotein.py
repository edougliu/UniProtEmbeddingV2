class UniProtein:
    def __init__(self, entry, entry_name, protein_name,
                 length, sequence, active_site, ecs_number, go_id):
        self.ecs_number = ecs_number
        self.entry_name = entry_name
        self.protein_name = protein_name
        self.gene_names = None
        self.organism = None
        self.length = length
        self.sequence = sequence
        self.active_site = active_site
        self.embedding = None  # there seem to be more unitpro proteins than embeddings
        self.entry = entry
        self.go_id = go_id

    def set_embedding(self, embedding):
        self.embedding = embedding
