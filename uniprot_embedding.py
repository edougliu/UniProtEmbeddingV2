import math

import numpy as np
import h5py
import random
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import sklearn
from sklearn.metrics.pairwise import cosine_similarity, euclidean_distances
import os
from scipy.spatial.distance import cityblock
from Bio import Align
from scipy.stats import pearsonr
# for later...
from goatools.base import download_go_basic_obo, download_ncbi_associations
from goatools.obo_parser import GODag
from goatools.go_enrichment import GOEnrichmentStudy
from goatools.anno.genetogo_reader import Gene2GoReader


# The raw embeddings (a .h5 file) are accessed here: https://www.uniprot.org/help/downloads#embeddings.
# The protein data ( a .tsv file) identified by UniProt ID are here:
# https://www.uniprot.org/uniprotkb?facets=reviewed%3Atrue&query=%2A.
# Click "download" then choose .tsv file and
# select what columns are needed

# Some Ideas:
#  >use Search Peptides analysis tool to select a set of proteins with an exact peptide match
#       and compare distances -> could a given parameter in the embeddings be mapped to a common peptide?


# Readme: UniProt is providing raw embeddings for UniProtKB/Swiss-Prot
# and some reference proteomes of model organisms.
# The embeddings are generated using the ProtT5 protein language model
# and stored in the standard HDF5 file format.
# There are two embeddings files generated: per-protein embeddings,
# where a fixed-length embeddings vector is computed for the whole protein sequence,
# and per-residue embeddings where a fixed-length embeddings vector is computed for each single residue.
# Note: Protein sequences longer than 12k residues are excluded due to limitation of GPU memory
# (this concerns only a handful of proteins).
def read_ecs(filename):
    ec_d = {}  # dictionary
    with open(filename) as f:
        header = next(f)
        for line in f:
            ws = line.split('\t')
            # element 0 is the protein accession number (id) and element 7 is the ecs number
            ec_d[ws[0]] = ws[7].strip().split(";")
    return ec_d


def read_active_site(filename):
    as_d = {}  # dictionary
    with open(filename) as f:
        header = next(f)
        for line in f:
            ws = line.split('\t')
            # element 0 is the protein accession number (id) and element 8 is the active site number
            as_d[ws[0]] = ws[7].strip().split(";")
    return as_d


def read_sequence(filename):
    seq_d = {}  # dictionary
    with open(filename) as f:
        header = next(f)
        for line in f:
            ws = line.split('\t')
            # element 0 is the protein id (protein accession) and element 7 is the sequence
            seq_d[ws[0]] = ws[7]
    return seq_d


def read_blast(filename):
    blast_d = {}  # dictionary
    with open(filename) as f:
        header = next(f)
        for line in f:
            ws = line.split('\t')
            # element 2 is the protein id (protein accession) and element 6 is the alignment score
            blast_d[ws[2]] = ws[6]
    return blast_d


def read_protein_length(filename, ec_num):
    ec_d = {}  # dictionary
    max_length = 0  # assuming length in the file refers to number of amino acids
    min_length = 34350  # num of amino acids in titin, largest protein
    total_length = 0
    avg_length = 0
    with open(filename) as f:
        header = next(f)
        for line in f:
            ws = line.split('\t')
            if ec_num in ws[7]:  # contains?
                length = int(ws[6])
                total_length = length + total_length
                ec_d[ws[0]] = length
                # print(ws[0] + ":" + str(length))
                if length > max_length:
                    max_length = length
                if length < min_length:
                    min_length = length
        avg_length = total_length / len(ec_d)
        # print("min length = " + str(min_length))
        # print("max length = " + str(max_length))
        # print("avg length = " + str(avg_length))
    return ec_d


def read_embeddings(pdict, proportion):
    # UniProt is providing raw embeddings (per-protein and per-residue using the ProtT5 model)
    # for UniProtKB/Swiss-Prot and some reference proteomes of model organisms
    # (such as Homo sapiens, SARS-CoV-2, and E. coli).
    # You can retrieve them from the Downloads page.
    with h5py.File("per-protein.h5", "r") as file:
        # print(f"number of entries: {len(file.items())}")
        for sequence_id, embedding in file.items():
            i = random.random()
            if i < proportion:
                pdict[sequence_id] = np.array(embedding)


def read_embeddings_blast(pdict, blast_dict):
    # load up the "embeddings" (vector containing all the parameters/features as a result of training)
    # associated with the proteins based on a set of proteins deemed very similar
    # via BLAST multi-sequence alignment score into a dictionary with key=protein id
    # and the value=the embedding
    # most of the proteins from the blast alignment file have no available raw embeddings
    # file2 = 'per-protein-subset-blast.h5'
    with h5py.File("per-protein.h5", "r") as file:
        for sequence_id in blast_dict:
            # print("looking for sequence id: " + sequence_id)
            try:
                embedding = file[sequence_id]
                # print("sequence found: " + sequence_id + " with embedding:" + str(embedding))
                pdict[sequence_id] = np.array(embedding)
            except KeyError:
                pass


def read_embeddings_ec(pdict, ec_num, ec_dict):
    # load up the "embeddings" (vector containing all the parameters/features as a result of training)
    # associated with the proteins for the specified enzyme category into a dictionary with key=protein id
    # and the value=the embedding
    file2 = 'per-protein-subset.h5'
    with h5py.File("per-protein.h5", "r") as file:
        # print(f"number of entries: {len(file.items())}")
        for sequence_id, embedding in file.items():
            if sequence_id in ec_dict:
                for num in ec_dict[sequence_id]:  # multiple ec nums
                    if ec_num in num:
                        # if num.startswith(ec_num):
                        create_subset_h5(file, file2, sequence_id)
                        # print("embedding:" + str(embedding))
                        pdict[sequence_id] = np.array(embedding)
                        # got one, so ignore others for this protein
                        break
            else:
                sys.stderr.write("not found:" + sequence_id + "\n")


def read_embeddings_ec2(pdict, ec_num1, ec_num2, ec_dict):
    with h5py.File("per-protein.h5", "r") as file:
        # print(f"number of entries: {len(file.items())}")
        for sequence_id, embedding in file.items():
            if sequence_id in ec_dict:
                for num in ec_dict[sequence_id]:  # multiple ec nums
                    if num.startswith(ec_num1) or num.startswith(ec_num2):
                        pdict[sequence_id] = np.array(embedding)
                        # got one, so ignore others for this protein
                        break
            else:
                sys.stderr.write("not found:" + sequence_id + "\n")


def calculate_embeddings_using_ec_num(ec_num_1, ec_num_2=None):
    # TODO refactor such that we can compare on any enzyme categorization available from UniProt
    #  (eg active site) Actually most of the enzyme categorizations are not given in a codified way
    #   that using them in a comparison is very difficult.
    #   Try using active site by extracting the amino acid at the active site and using that to compare
    #   Do proteins with same amino acid at active site have similar embeddings?
    # Read in the file that maps protein accession to EC numbers
    # Can have multiple EC numbers per protein.
    ec_dict = read_ecs(
        "uniprot-compressed_true_download_true_fields_accession_2Creviewed_2C-2023.03.08-23.38.08.49.tsv")

    # this is the file with the additional attributes
    # as_dict = read_active_site(
    #   "uniprot-compressed_true_download_true_fields_accession_2Cid_2Cprotei-2023.06.26-03.19.53.58.tsv")
    seq_dict = read_sequence(
        "uniprotkb_with_sequences_2023_07_03.tsv")

    # What proportion of the data to use (e.g. 0.01 = 1%)
    proportion = 0.01

    # Mapping from a protein accession to embedding data
    pdict = {}
    # Mapping from protein accession to sequence

    # Read in suitable embeddings (for given ec numbers or for a random subset)
    if ec_num_2 is None:
        read_embeddings_ec(pdict, ec_num_1, ec_dict)
        file_prefix = ec_num_1
    else:
        read_embeddings_ec2(pdict, ec_num_1, ec_num_2, ec_dict)
        file_prefix = ec_num_1 + ec_num_2
    # or
    # read_embeddings(pdict, proportion) # uses just a proportion
    read_embeddings(pdict, file_prefix, seq_dict)


def calculate_embeddings_using_active_site(sdict, as_dict, ids, i, j):
    seq = sdict[ids[j]]
    active_site = as_dict[ids[j]]
    print(active_site)
    print(len(active_site))
    if active_site[0] != '':
        active_site_num_str = active_site[0]  # strip off ACT_SITE, there can be more than one of these
        active_site_num = active_site_num_str[9:]
        print("active site: <" + active_site[0] + ">")
        print("active site num:" + active_site_num)
        print(len(active_site_num))
        active_site_aa = seq[int(active_site_num)]
        print("active site amino acid: " + active_site_aa)
        return active_site_aa
    else:
        return None


def calculate_embeddings_using_blast():
    # used https://www.uniprot.org/blast
    # took the first protein from the original tsv file (A0A023I7E1)
    # and used that to get the top matches via sequence alignment (ncbiblast-R20230626-032626-0211-47002547-p1m.tsv).
    # use these to calculate embeddings to see how the embedding distances from a set of sequences known to be similar
    # per BLAST multi-sequence alignment compare
    blast_dict = read_blast("ncbiblast-R20230626-032626-0211-47002547-p1m.tsv")  # top similar sequences per blast
    # Mapping from a protein accession to embedding data
    # check to see if you can select just SwissProt proteins
    pdict = {}
    read_embeddings_blast(pdict, blast_dict)
    file_prefix = "blast"
    read_embeddings(pdict, file_prefix)  # not sure why I have this-has wrong num of parameters


def read_embeddings(pdict, file_prefix, sdict, as_dict=None):
    # How many did we read in?
    ids = list(pdict.keys())
    # the number of proteins that we found associated with the given enzyme category
    num_proteins = len(ids)
    print("num of proteins found " + str(num_proteins))

    # lets make a large array for the all-vs-all distances
    # but I was being cheap and not doing a 2D array, because we'd only
    # need to fill half of it and the other half would be the same
    # as d(x,y) = d(y,x) so just doing a 1D array

    # compare the protein embeddings for all the proteins associated with the given enzyme category
    # to determine the similarity as represented by the euclidean distance
    #  variable "distances" will store the distance between each protein and each of the other proteins
    #  after it for graphing
    distance_norm = np.empty(int(num_proteins * (num_proteins - 1) / 2))  # ndarray
    distance_cosine = np.empty(int(num_proteins * (num_proteins - 1) / 2))  # ndarray
    distance_euclidean = np.empty(int(num_proteins * (num_proteins - 1) / 2))  # ndarray
    distance_manhattan = np.empty(int(num_proteins * (num_proteins - 1) / 2))  # ndarray
    distance_alignment = np.empty(int(num_proteins * (num_proteins - 1) / 2))  # ndarray
    distance_pearson = np.empty(int(num_proteins * (num_proteins - 1) / 2))  # ndarray

    max_norm = 0
    max_euclidean = 0
    max_manhattan = 0
    max_cosine = 0
    max_alignment = 0
    max_pearson = 0

    k = 0  # next position in array to fill
    for i in range(num_proteins):
        for j in range(i + 1, num_proteins):
            # print(i,j) # we'll calculate distance between proteins i and j the norm of a matrix is a measure of the
            # "size" that is unrelated to its number of rows and columns in linear algebra, the euclidean norm is the
            # square root of the sum of all the squares so d is a similarity score indicating how similar the given
            # protein embedding is to all the others

            # The cosine similarity always belongs to the interval [-1,
            # 1] For example, two proportional vectors have a cosine similarity of 1, two orthogonal vectors have a
            # similarity of 0, and two opposite vectors have a similarity of -1. In some contexts, the component
            # values of the vectors cannot be negative, in which case the cosine similarity is bounded in [0,
            # 1] For example, in information retrieval and text mining, each word is assigned a different coordinate
            # and a document is represented by the vector of the numbers of occurrences of each word in the document.
            # Cosine similarity then gives a useful measure of how similar two documents are likely to be,
            # in terms of their subject matter, and independently of the length of the documents.
            # cosine has the nice property that it is 1.0 for identical vectors and 0.0 for orthogonal vectors)

            # print("original shape: " + str(pdict[ids[i]].shape))
            # print("calculating distance of " + str(pdict[ids[i]].reshape(1, -1)))
            # print(" and " + str(pdict[ids[j]].reshape(1, -1)))

            # TODO build on this to be able to select a set of proteins with same amino acid at active site
            # active_site_aa = calculate_embeddings_using_active_site(sdict, as_dict, ids, i, j)

            # this is the same as euclidean
            # method_of_comparison = np.linalg.norm
            distance_n = calculate_norm(pdict, ids, i, j)

            # method_of_comparison = cosine_similarity
            distance_c = calculate_cosine(pdict, ids, i, j)

            # method_of_comparison = euclidean_distances
            distance_e = calculate_euclidean(pdict, ids, i, j)

            # method_of_comparison = cityblock  # manhattan
            distance_m = calculate_manhattan(pdict, ids, i, j)  # x values

            # method_of_comparison = alignment  # a high score means a more similar sequence
            distance_a = calculate_alignment(sdict, ids, i, j)  # y value

            # method_of_comparison = Pearson correlation
            distance_p = calculate_pearson(pdict, ids, i,
                                           j)  # <- this should be performed on the above as two lists of values

            distance_norm[k] = distance_n
            if distance_norm[k] > max_norm:
                max_norm = distance_norm[k]
                # print("max_norm is: " + str(max_norm))

            distance_cosine[k] = distance_c
            if distance_cosine[k] > max_cosine:
                max_cosine = distance_cosine[k]
                # print("max_cosine is: " + str(max_cosine))

            distance_euclidean[k] = distance_e
            if distance_cosine[k] > max_euclidean:
                max_euclidean = distance_euclidean[k]
                # print("max_euclidean is: " + str(max_euclidean))

            distance_manhattan[k] = distance_m
            if distance_manhattan[k] > max_manhattan:
                max_manhattan = distance_manhattan[k]
                # print("max_manhattan is: " + str(max_manhattan))

            distance_alignment[k] = distance_a
            if distance_alignment[k] > max_alignment:
                max_alignment = distance_alignment[k]
                # print("max_alignment is: " + str(max_alignment))

            distance_pearson[k] = distance_p
            if distance_alignment[k] > max_pearson:
                max_pearson = distance_pearson[k]
                # print("max_pearson is: " + str(max_pearson))
            k += 1
    print(distance_euclidean)
    plot(distance_norm, distance_alignment, "norm", file_prefix, math.ceil(max_norm))
    # cosine should be the range between 0 and 1
    # but for some reason there are a small number of datapoints past 1
    plot(distance_cosine, distance_alignment,  "cosine", file_prefix, math.ceil(max_cosine))
    plot(distance_euclidean, distance_alignment, "euclidean", file_prefix, math.ceil(max_euclidean))
    plot(distance_manhattan, distance_alignment, "manhattan", file_prefix, math.ceil(max_manhattan))
    # plot(distance_alignment, "alignment", file_prefix, math.ceil(max_alignment), 800)
    plot(distance_pearson, distance_alignment, "pearson", file_prefix, math.ceil(max_pearson))


def plot(distances, alignment_scores, method_of_comparison, ec_num, max_range=6, max_domain=None):
    # Plot embedding data
    #  kde = kernal density estimates,
    #  a non-parametric (no underlying assumptions about shape of data, eg bell shaped)
    #  statistical technique used to estimate the probability density function (PDF) of a random variable
    # sns.histplot(distances, binwidth=0.25, kde=True) # change this to a scatter plot

    # change to plot blast similarity scores with the embedding distance
    # indices = np.triu_indices(pairwise_distances.shape[0], k=1)  we already have the lower triangular indices
    # indices = np.indices(distances.shape) we have just a 1D array
    # indices = np.arange(len(distances))
    plt.figure(figsize=(8, 6))
    plt.scatter(distances, alignment_scores)
    plt.xlabel('pair-wise distance')
    plt.ylabel('blast alignment score')
    plt.title('Pairwise Embedding Distances Scatter Plot using ' + method_of_comparison)
    plt.savefig("ec_" + ec_num + "_" + method_of_comparison + ".png", dpi=400)
    plt.close()


def calculate_pearson(pdict, ids, i, j):
    p = pearsonr(pdict[ids[i]], pdict[ids[j]])
    # print(p.statistic)
    return p.statistic


def calculate_alignment(sdict, ids, i, j):
    sequence1 = sdict[ids[i]]  # look up the sequence for the given protein id
    sequence2 = sdict[ids[j]]
    a = perform_sequence_alignment(sequence1, sequence2)
    return a


def calculate_manhattan(pdict, ids, i, j):
    # Σ|Ai – Bi|
    d = cityblock(pdict[ids[i]], pdict[ids[j]])
    return d


def calculate_euclidean(pdict, ids, i, j):
    # d = method_of_comparison(pdict[ids[i]] - pdict[ids[j]])  # euclidean ?
    # reshape reduces the number of dimensions of a matrix
    # (Σ(Ai-Bi)^2)^1/2
    d = euclidean_distances(pdict[ids[i]].reshape(1, -1), pdict[ids[j]].reshape(1, -1))
    return d


# a vector norm assigns a non-negative value to a vector, representing its length
# or magnitude in a vector space. Different norms, such as the Euclidean norm, Manhattan norm,
# and maximum norm, capture different notions of "size" or "distance"
def calculate_norm(pdict, ids, i, j):
    d = np.linalg.norm(pdict[ids[i]] - pdict[ids[j]])
    # d = np.linalg.norm(pdict[ids[i]].reshape(1, -1), pdict[ids[j]].reshape(1, -1))
    return d


def calculate_cosine(pdict, ids, i, j):
    d = cosine_similarity(pdict[ids[i]].reshape(1, -1), pdict[ids[j]].reshape(1, -1))
    # d = cosine_similarity(pdict[ids[i]], pdict[ids[j]])
    return d


def calculate_simple_length(ec_num):
    # TODO de-dup
    protein_length_dict = read_protein_length(
        "uniprot-compressed_true_download_true_fields_accession_2Creviewed_2C-2023.03.08-23.38.08.49.tsv", ec_num)

    # What proportion of the data to use (e.g. 0.01 = 1%)
    proportion = 0.01

    # Mapping from a protein accession to embedding data
    pdict = protein_length_dict

    # How many did we read in?
    ids = list(pdict.keys())
    # the number of proteins that we found associated with the given enzyme category
    num_proteins = len(ids)
    # print("num of proteins found for enzyme category " + ec_num + ": " + str(num_proteins))
    distances = np.empty(int(num_proteins * (num_proteins - 1)))  # ndarray
    lengths = np.empty(int(num_proteins))
    k = 0  # next position in array to fill
    for i in range(num_proteins):
        lengths[i] = pdict[ids[i]]
        # print("length:" + str(lengths[i]))
        for j in range(i + 1, num_proteins):
            d = abs(pdict[ids[i]] - pdict[ids[j]])  # simple diff between the length of the proteins
            # print("distance: " + str(d))
            distances[k] = abs(d)
            k += 1

    # plot lengths
    sns.displot(lengths, binwidth=0.25, kde=True)
    plt.xlim(0, 6)
    plt.xticks(range(6))
    plt.savefig("ec_protein_length_" + ec_num + ".png", dpi=400)
    plt.close()

    # Plot differences in lengths data
    sns.displot(distances, binwidth=0.25, kde=True)
    plt.xlim(0, 6)  # why is this 6?
    plt.xticks(range(6))
    plt.xlabel('pair-wise distance')
    plt.ylabel('data points')
    plt.title('Histogram of embedding distances')
    plt.savefig("ec_protein_length_diff_" + ec_num + ".png", dpi=400)
    plt.close()
    # or
    # plt.savefig("ec_"+ec_num1+"_"+ec_num2+".png",dpi=400)
    # or
    # plt.savefig("random"+f'{proportion:.3f}'+"distn_3.png",dpi=400)
    # Plot simple distance data
    # d = method_of_comparison(pdict[ids[i]], pdict[ids[j]])


# trying to stream the subset of protein data that is selected
# for analysis to a separate file so as to be inspect-able using HDFViewer
def create_subset_h5(larger_h5, subset_h5, dataset_name):
    small_h5 = h5py.File(subset_h5, 'a')  # 'per-protein-subset.h5'
    larger_h5.copy(dataset_name, small_h5)
    # print("dataset name=" + dataset_name)


def perform_sequence_alignment(sequence1, sequence2):
    # BiopythonDeprecationWarning: Bio.pairwise2 has been deprecated, and we intend to remove it in a future release
    # of Biopython. As an alternative, please consider using Bio.Align.PairwiseAligner as a replacement, and contact
    # the Biopython developers if you still need the Bio.pairwise2 module.
    # alignments = pairwise2.align.globalxx(sequence1, sequence2)
    # alignment_score = alignments[0].score
    # print("Alignment Score:", alignment_score)
    # return alignment_score
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
    alignments = aligner.align(sequence1, sequence2)
    # print("Alignment Score:", alignments.score)
    return alignments.score


def main():
    # TODO try comparing within other categories such as active site (need to download)
    # try comparing to a completely different enzyme category
    ec_num = "5.1.1.1"
    ec_num_1 = "4.1.1.1"
    ec_num_2 = "5.1.1.1"
    # before running clean up old subset h5 file
    file2 = 'per-protein-subset.h5'
    if os.path.exists(file2):
        os.remove(file2)
    # uncomment one or more of the below options:
    calculate_embeddings_using_ec_num(ec_num)
    # calculate_embeddings_using_ec_num(ec_num_1, ec_num_2)
    # calculate_embeddings_using_blast()
    # calculate_simple_length(ec_num_1)
    # calculate_embeddings_using_active_site()


main()
