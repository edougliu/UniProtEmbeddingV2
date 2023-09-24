import os
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics.pairwise import cosine_similarity, euclidean_distances
from Bio import Align
from scipy.stats import pearsonr
from Bio.Align import substitution_matrices
import subprocess

from diamond_result import DiamondResult


class Analyser:
    def analyse(self, uniproteins, file_prefix):
        print("Analyser current working directory:", os.getcwd())
        num_proteins = len(uniproteins)
        print("num of proteins found " + str(num_proteins))

        distances_euclidean = np.empty(int(num_proteins * (num_proteins - 1) / 2))  # ndarray
        alignments = np.empty(int(num_proteins * (num_proteins - 1) / 2))  # ndarray
        diamonds = np.empty(int(num_proteins * (num_proteins - 1) / 2))  # ndarray

        k = 0
        for i in range(num_proteins):
            for j in range(i + 1, num_proteins):
                uniprotein_curr = uniproteins[i]
                uniprotein_next = uniproteins[j]
                distance = self.calculate_euclidean(uniprotein_curr, uniprotein_next)
                alignment = self.calculate_alignment(uniprotein_curr, uniprotein_next)
                diamond = self.calculate_diamond(uniprotein_curr, uniprotein_next)
                distances_euclidean[k] = distance
                alignments[k] = alignment
                diamonds[k] = diamond
                k += 1
        pearson_correlation = self.calculate_pearson(distances_euclidean, alignments)
        note = "pearson correlation: " + str(pearson_correlation)
        print(note)
        self.plot(distances_euclidean, alignments, "euclidean", file_prefix, note)

    def plot(self, distances, alignment_scores, method_of_comparison, selection_criteria, note):
        plt.figure(figsize=(8, 6))
        plt.scatter(distances, alignment_scores)
        plt.text(
            0.5,  # x-coordinate of the text
            1.5,  # y-coordinate of the text (adjusted to be below the title)
            note,  # got tired of fiddling with this to try to make it show up on the plot
            fontsize=12,
            color='red',
            ha='center',  # Horizontal alignment ('center' aligns the text with the plot title)
        )
        plt.xlabel('pair-wise distance')
        plt.ylabel('alignment score')
        plt.title('Pairwise Embedding Distances Scatter Plot using ' + method_of_comparison)
        plt.savefig("distances_for_" + selection_criteria + "_" + method_of_comparison + ".png", dpi=400)
        plt.close()

    def calculate_pearson(self, distances, alignments):
        print("calculating pearson")
        print(distances)
        print(alignments)
        p = pearsonr(distances, alignments)
        return p.statistic

    def calculate_alignment(self, uniprotein_curr, uniprotein_next):
        sequence1 = uniprotein_curr.sequence
        sequence2 = uniprotein_next.sequence
        a = self.perform_sequence_alignment(sequence1, sequence2)
        return a

    def calculate_diamond(self, uniprotein_curr, uniprotein_next):
        sequence1 = uniprotein_curr.sequence
        sequence2 = uniprotein_next.sequence
        d = self.perform_diamond_similarity(sequence1, sequence2)
        # turn the list into a nice object
        accession = d[0]
        diamond = DiamondResult.create_instance_from_list(d)
        print("comparing " + uniprotein_curr.entry + " with "
              + uniprotein_next.entry + " bit score:" + str(diamond.bit_score))
        return diamond.bit_score

    def calculate_euclidean(self, uniprotein_curr, uniprotein_next):
        # reshape reduces the number of dimensions of a matrix
        # (Σ(Ai-Bi)^2)^1/2
        # euclidean distance of 0 means identical and the higher the score, the more different
        d = euclidean_distances(uniprotein_curr.embedding.reshape(1, -1),
                                uniprotein_next.embedding.reshape(1, -1))
        return d

    def perform_sequence_alignment(self, sequence1, sequence2):
        # set mode to use local and report the % that aligned as well as the score
        # look at Bit score from Diamond (a good match would be ~60)
        # try removing outliers so focus in on the "clump"
        # are alignment scores are based on length? normalize for length but Diamond is the better option
        # Diamond: take all sequences, create database from them, use other diamond command to search
        # each one again the database you build,
        # outputs file, each line in the file is a pair of sequences' alignment score
        # https: //github.com/bbuchfink/diamond/wiki
        # download sequences in fasta format instead of tsv
        aligner = Align.PairwiseAligner()
        blosum62_matrix = substitution_matrices.load("BLOSUM62")
        aligner.substitution_matrix = blosum62_matrix
        aligner.mode = "local"
        #  no idea what to set these parameters to or what they might mean...
        aligner.match_score = 2
        aligner.mismatch_score = -1
        aligner.open_gap_score = -5
        aligner.extend_gap_score = -2
        alignments = aligner.align(sequence1, sequence2)
        return alignments.score



    def perform_diamond_similarity(self, sequence1, sequence2):
        # use a previously created diamond database and the diamond executable
        # to perform a comparison between two sequences
        """
        From the documentation at https://github.com/bbuchfink/diamond/wiki/1.-Tutorial
        The file is generated in tabular-separated (TSV) format composed of 12 fields, a format corresponding to the
        format generated by BLAST using the option -outfmt 6. The 12 fields are:
        ->Query accession: the accession of
        the sequence that was the search query against the database, as specified in the input FASTA file after the >
        character until the first blank.
        ->Target accession: the accession of the target database sequence (also called
        subject) that the query was aligned against.
        ->Sequence identity: The percentage of identical amino acid
        residues that were aligned against each other in the local alignment.
        ->Length: The total length of the local alignment, which including matching and mismatching positions of query
        and subject, as well as gap positions in the query and subject.
        ->Mismatches: The number of non-identical amino acid residues aligned against each other.
        ->Gap openings: The number of gap openings.
        ->Query start: The starting coordinate of the local alignment in the query (1-based).
        ->Query end: The ending coordinate of the local alignment in the query (1-based).
        ->Target start: The starting coordinate of the local alignment in the target (1-based).
        ->Target end: The ending coordinate of the local alignment in the target (1-based).
        ->E-value: The expected value of the hit quantifies
        the number of alignments of similar or better quality that you expect to find searching this query against a
        database of random sequences the same size as the actual target database. This number is most useful for
        measuring the significance of a hit. By default, DIAMOND will report all alignments with e-value < 0.001,
        meaning that a hit of this quality will be found by chance on average once per 1,000 queries.
        ->Bit score: The bit score is a scoring matrix independent measure of the (local) similarity
        of the two aligned sequences, with higher numbers meaning more similar.
        It is always >= 0 for local Smith Waterman alignments.
        """
        diamond_executable_path = "C:\\myProjects\\Tools\\diamond-windows\\diamond.exe"
        diamond_database_path = "C:\\myProjects\\Tools\\uniprot_sprot.fasta\\uniprot.dmnd"
        diamond_command = [
            diamond_executable_path,
            "blastp",  # Use "blastp" for protein-protein comparison
            "-d", diamond_database_path,
            "--quiet",  # Suppress DIAMOND's progress output
            "--query", "-",  # Read the query sequence from stdin
            "--outfmt", "6",  # Specify the output format (tabular format)
        ]

        with subprocess.Popen(
                diamond_command,
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True,  # Use text mode for communication
                shell=True
        ) as process:
            # Write the sequences to DIAMOND's stdin
            query_data = f">sequence1\n{sequence1}\n>sequence2\n{sequence2}\n"
            stdout, stderr = process.communicate(input=query_data)
        print("stdout: " + stdout)
        if stderr:
            print("DIAMOND Error:", stderr)

        # stdout is just a string tokenized with tabs so turn it into a list
        return stdout.split('\t')
