import numpy as np
import pandas as pd
import h5py
import random
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import sklearn
from sklearn.metrics.pairwise import cosine_similarity, euclidean_distances
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA


def biplot(title, score, coef, labels=None):
    plt.title(title)
    xs = score[:, 0]
    ys = score[:, 1]
    n = coef.shape[0]
    scalex = 1.0 / (xs.max() - xs.min())
    scaley = 1.0 / (ys.max() - ys.min())
    plt.scatter(xs * scalex, ys * scaley,
                s=5,
                color='orange')

    for i in range(n):
        plt.arrow(0, 0, coef[i, 0],
                  coef[i, 1], color='purple',
                  alpha=0.5)
        plt.text(coef[i, 0] * 1.15,
                 coef[i, 1] * 1.15,
                 labels[i],  # <-error here - make up some labels?
                 color='darkblue',
                 ha='center',
                 va='center')

    plt.xlabel("PC{}".format(1))
    plt.ylabel("PC{}".format(2))

    # plt.figure()
    plt.show()


# number of entries: 568354
# embedding size: 1024
def read_ecs(filename):
    ec_d = {}  # dictionary
    with open(filename) as f:
        header = next(f)
        for line in f:
            ws = line.split('\t')
            ec_d[ws[0]] = ws[7].strip().split(";")
    return ec_d


def read_embeddings(pdict, proportion):
    # UniProt is providing raw embeddings (per-protein and per-residue using the ProtT5 model)
    # for UniProtKB/Swiss-Prot and some reference proteomes of model organisms
    # (such as Homo sapiens, SARS-CoV-2, and E. coli).
    # You can retrieve them from our Downloads page.
    with h5py.File("per-protein.h5", "r") as file:
        # print(f"number of entries: {len(file.items())}")
        for sequence_id, embedding in file.items():
            i = random.random()
            if i < proportion:
                pdict[sequence_id] = np.array(embedding)


def read_embeddings_ec(pdict, ec_num, ec_dict):
    # load up the "embeddings" (vector containing all the parameters/features as a result of training)
    # associated with the proteins for the specified enzyme category into a dictionary with key=protein id
    # and the value=the embedding
    with h5py.File("per-protein.h5", "r") as file:
        # print(f"number of entries: {len(file.items())}")
        for sequence_id, embedding in file.items():
            if sequence_id in ec_dict:
                for num in ec_dict[sequence_id]:  # multiple ec nums
                    if ec_num in num:
                        # if num.startswith(ec_num):
                        pdict[sequence_id] = np.array(embedding)
                        # got one, so ignore others for this protein
                        break
            else:
                sys.stderr.write("not found:" + sequence_id + "\n")


def perform_pca(ec_num):
    # Read in the file that maps protein accession to EC numbers
    # Can have multiple EC numbers per protein.
    ec_dict = read_ecs(
        "uniprot-compressed_true_download_true_fields_accession_2Creviewed_2C-2023.03.08-23.38.08.49.tsv")

    # What proportion of the data to use (e.g. 0.01 = 1%)
    proportion = 0.01

    # Mapping from a protein accession to embedding data
    pdict = {}

    # Read in suitable embeddings (for given ec numbers or for a random subset)
    read_embeddings_ec(pdict, ec_num, ec_dict)  # uses all, not a proportion

    # How many did we read in?
    ids = list(pdict.keys())
    # the number of proteins that we found associated with the given enzyme category
    num_proteins = len(ids)
    print("num of proteins found for enzyme category " + ec_num + ": " + str(num_proteins))

    df = pd.DataFrame(data=pdict.values(), columns=range(1, pdict.__sizeof__()))
    df.head(6)
    print("the embedding dataframe has size:" + str(df.size))
    scaler = StandardScaler()
    scaler.fit(df)
    protein_embeddings_scaled = scaler.transform(df)
    dataframe_scaled = pd.DataFrame(data=protein_embeddings_scaled)
    dataframe_scaled.head(6)
    print("the scaled embedding dataframe has size:" + str(dataframe_scaled.size))
    # let's just try 10 to start with TODO: optimize
    pca = PCA(n_components=10)
    pca.fit_transform(dataframe_scaled)
    prop_var = pca.explained_variance_ratio_
    eigenvalues = pca.explained_variance_
    pc_numbers = np.arange(pca.n_components_) + 1
    # A scree plot is a graphical tool used in exploratory factor analysis (EFA)
    # or principal component analysis (PCA) to assess the number of meaningful factors
    # or components in a dataset. It helps determine the optimal number of factors
    # or components to retain.
    # ie What is the "right" number of components to retain where by you can retain most of the important information
    plt.plot(pc_numbers,
             prop_var,
             'ro-')
    plt.title('Figure 1: Scree Plot', fontsize=8)
    plt.ylabel('Proportion of Variance', fontsize=8)
    plt.savefig("ec_protein_pca_" + ec_num + ".png")
    plt.close()
    # re-run PCA with first two components per the elbow rule
    pca = PCA(n_components=2)
    PC = pca.fit_transform(dataframe_scaled)
    pca_low_dim = pd.DataFrame(data=PC,
                               columns=['PC1', 'PC2'])
    pca_low_dim.head(6)  # this doesn't work
    title = 'Biplot of PCA for Proteins for enzymes ' + ec_num
    biplot(title, PC,
           np.transpose(pca.components_),
           list(pdict.feature_names))  # <-we don't know what our features are


def main():
    # taken from tutorial at https://statisticsglobe.com/principal-component-analysis-python
    ec_num = "5.1.1.1"
    # ec_num1 = "4.1.1.1"
    # ec_num2 = "5.1.1.1"
    perform_pca(ec_num)


main()
