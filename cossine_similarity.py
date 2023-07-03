import numpy as np
import funcoes_gerais as f
import scipy.cluster.hierarchy as shc
import matplotlib.pyplot as plt
from numpy.linalg import norm
from scipy.cluster.hierarchy import single, cophenet
from scipy.spatial.distance import pdist, squareform

def cosine_similarity(vector1, vector2):
    """
    Given two vectors, return the cosine similarity between them.
    """
    vector1 = np.array(vector1)
    vector2 = np.array(vector2)
    similarity = np.dot(vector1,vector2)/(norm(vector1)*norm(vector2))
    return similarity

def word_based_method_cos(seq1,seq2,k):
    """
    Given two sequence-representing vectors and a k, returns the Cosine Similarity between them..
    """
    
    w1 = f.generate_kmers(seq1, k)
    w2 = f.generate_kmers(seq2, k)
    
    w_union = f.union_vectors(w1, w2)
    
    count1= f.count_items(w_union, w1)
    count2= f.count_items(w_union, w2)

    distance = cosine_similarity(count1,count2)
    
    return distance


def best_k_cos(seq1,seq2):
    
    """
    Given two sequence-representing vectors, returns the k that minimizes the Cosine Similarity between them (and the Cosine Similarity)
    """   
    
    kmax=15
    dist = word_based_method_cos(seq1,seq2,1)
    k=1
    
    for i in range(2,kmax + 1):
        if dist >word_based_method_cos(seq1,seq2,i):
            k = i
            dist = word_based_method_cos(seq1,seq2,i)

    return k,dist



def all_k_cos(seq1,seq2):
    kmax=15
    dist = []
    k =list(range(1,kmax+1))
    
    for i in range(1,kmax+1):
        dist += [word_based_method_cos(seq1,seq2,i)]
    
    return k,dist



#FOR SEVERAL SEQUENCES (AND NOT JUST 2)

def matrix_k_cos(vector):
    
    """
    Given a vector with n sequences,creates a matrix where M(a,b) is the best k for the cossine similarity
    """
    
    n =len(vector)
    M = np.zeros((n,n))
    for i in range(n):
        for j in range(i,n):
            M[i,j] = best_k_cos(vector[i],vector[j])[0]
            
    median = int(np.median(M))
            
    return [median,M]

def matrix_dist_cos(vector):
    k = matrix_k_cos(vector)[0]
    n =len(vector)
    M = np.zeros((n,n))
    
    for i in range(n):
        for j in range(n):
            M[i,j] = word_based_method_cos(vector[i],vector[j],k)
            
    return M

def plot_cos(genes,labels, cluster):
    n_plots = len(cluster)
    
    distance_matrix = matrix_dist_cos(genes)
    
    if n_plots <=3:
        for i in range(n_plots):
            # calculate the linkage matrix
            linkage_matrix = shc.linkage(distance_matrix, method =  cluster[i])
            print(linkage_matrix)
            shc.set_link_color_palette(['#d81159' ,  '#8f2d56' ,   '#218380' ,   '#fbb13c' ,   '#73d2de'])

            # COPHENETIC COEFFICIENT
            Y = pdist(distance_matrix)
            Z = linkage_matrix
            c, coph_dists = cophenet(Z, Y)
            print(f"Cophenetic correlation coefficient for {cluster[i]}: {c}")

            # plot the dendrogram with horizontal orientation
            plt.subplot(1,n_plots,i+1)
            plt.title("Hierarquy Cluster:" + cluster[i], fontsize = 10, color = 'black')
            dendrogram = shc.dendrogram(linkage_matrix, orientation='left', labels =labels)
    
    if n_plots == 4:
        for i in range(n_plots):
            # calculate the linkage matrix
            linkage_matrix = shc.linkage(distance_matrix, method =  cluster[i])
            shc.set_link_color_palette(['#d81159' ,  '#8f2d56' ,   '#218380' ,   '#fbb13c' ,   '#73d2de'])

            # COPHENETIC COEFFICIENT
            Y = pdist(distance_matrix)
            Z = linkage_matrix
            c, coph_dists = cophenet(Z, Y)
            print(f"Cophenetic correlation coefficient for {cluster[i]}: {c}")

            # plot the dendrogram with horizontal orientation
            plt.subplot(2,2,i+1)
            plt.title("Hierarquy Cluster:" + cluster[i], fontsize = 10, color = 'black')
            dendrogram = shc.dendrogram(linkage_matrix, orientation='left', labels =labels)



    if n_plots > 4:
        for i in range(n_plots):
            # calculate the linkage matrix
            linkage_matrix = shc.linkage(distance_matrix, method =  cluster[i])
            shc.set_link_color_palette(['#d81159' ,  '#8f2d56' ,   '#218380' ,   '#fbb13c' ,   '#73d2de'])

            # COPHENETIC COEFFICIENT
            Y = pdist(distance_matrix)
            Z = linkage_matrix
            c, coph_dists = cophenet(Z, Y)
            print(f"Cophenetic correlation coefficient for {cluster[i]}: {c}")
            
            plt.subplot(2,3,i+1)
            plt.title("Hierarquy Cluster:  " + " " + cluster[i], fontsize = 10, color = 'black')
            dendrogram = shc.dendrogram(linkage_matrix, orientation='left', labels =labels)

    

    plt.suptitle("Cosine Similarity")
    plt.subplots_adjust(left=0.1, right=0.9, bottom=0.1, top=0.9, wspace=0.5, hspace=0.5)
    plt.subplots.figsize = (20,40)
    plt.show()




