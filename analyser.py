import os
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics.pairwise import cosine_similarity, euclidean_distances
from Bio import Align
from scipy.stats import pearsonr
from Bio.Align import substitution_matrices


class Analyser:
    def analyse(self, uniproteins, file_prefix):
        print("Analyser current working directory:", os.getcwd())
        num_proteins = len(uniproteins)
        print("num of proteins found " + str(num_proteins))

        distances_euclidean = np.empty(int(num_proteins * (num_proteins - 1) / 2))  # ndarray
        alignments = np.empty(int(num_proteins * (num_proteins - 1) / 2))  # ndarray

        k = 0
        for i in range(num_proteins):
            for j in range(i + 1, num_proteins):
                uniprotein_curr = uniproteins[i]
                uniprotein_next = uniproteins[j]
                distance = self.calculate_euclidean(uniprotein_curr, uniprotein_next)
                alignment = self.calculate_alignment(uniprotein_curr, uniprotein_next)
                distances_euclidean[k] = distance
                alignments[k] = alignment
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
            note,  # The text you want to add
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
