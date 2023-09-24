class DiamondResult:
    def __init__(self, accession, seq_identity_percent, length,
                 mismatches, gaps, query_start, query_end, target_start,
                 target_end, e_value, bit_score):
        self.accession = accession
        self.seq_identity_percent = seq_identity_percent
        self.length = length
        self.mismatches = mismatches
        self.gaps = gaps
        self.query_start = query_start
        self.query_end = query_end
        self.target_start = target_start
        self.target_end = target_end
        self.e_value = e_value
        self.bit_score = bit_score

    @classmethod  # because python overloading is lame
    def create_instance_from_list(cls, property_list):
        print("creating diamond result from " + str(property_list))
        return cls(property_list[1], property_list[2], property_list[3],
                   property_list[4], property_list[5], property_list[6],
                   property_list[7], property_list[8], property_list[9],
                   property_list[10], property_list[11])
