import numpy as np
import funcoes_gerais as f
import scipy.cluster.hierarchy as shc
import matplotlib.pyplot as plt
from math import sqrt
from scipy.cluster.hierarchy import single, cophenet
from scipy.spatial.distance import pdist, squareform

def manhattan(a, b):
    return sum(abs(val1-val2) for val1, val2 in zip(a,b))

# MANHATTAN DISTANCE -----------------------------------------------------------------------------------------------------------------------

def word_based_method_man(seq1,seq2,k):
    """
    Given two sequence-representing vectors and a k, returns the Manhattan distance between them.
    """
    
    w1 = f.generate_kmers(seq1, k)
    w2 = f.generate_kmers(seq2, k)
    
    w_union = f.union_vectors(w1, w2)
    
    count1= f.count_items(w_union, w1)
    count2= f.count_items(w_union, w2)

    distance = manhattan(count1,count2)
    
    return distance

def best_k_man(seq1,seq2):
    
    """
    Given two sequence-representing vectors, returns the k that minimizes the Manhattan distance between them (and the Manhattan distance itself)
    """   
    
    kmax=15 # ate 15
    dist = word_based_method_man(seq1,seq2,1)
    k=1
    
    for i in range(2,kmax + 1):
        if dist >word_based_method_man(seq1,seq2,i):
            k = i
            dist = word_based_method_man(seq1,seq2,i)

    return k,dist

def all_k_man(seq1,seq2):
    kmax=min(len(seq1),len(seq2))
    dist = []
    k =list(range(1,kmax+1))
    
    for i in range(1,kmax+1):
        dist += [word_based_method_man(seq1,seq2,i)]
    
    return k,dist


    return k,dist

# MANHATTAN DISTANCE FOR MORE THAN 2 SEQUENCES

def matrix_k_man(vector):
    
    """
    Given a vector with n sequences,creates a matrix where M(a,b) is the best k for the normalized Manhattan distance
    """
    
    n =len(vector)
    M = np.zeros((n,n))
    for i in range(n):
        for j in range(i,n):
            M[i,j] = best_k_man(vector[i],vector[j])[0]
            
    median = int(np.median(M))
            
    return [median,M]

def matrix_dist_man(vector):
    k = matrix_k_man(vector)[0]
    n =len(vector)
    M = np.zeros((n,n))
    
    for i in range(n):
        for j in range(n):
            M[i,j] = word_based_method_man(vector[i],vector[j],k)
            
    return M

def plot_man(genes,labels,cluster):
    n_plots = len(cluster)
    distance_matrix_man = matrix_dist_man(genes)

    if n_plots <=3:
        for i in range(n_plots):
            # calculate the linkage matrix
            linkage_matrix = shc.linkage(distance_matrix_man, method =  cluster[i])
            shc.set_link_color_palette(['#d81159' ,  '#8f2d56' ,   '#218380' ,   '#fbb13c' ,   '#73d2de'])

            # COPHENETIC COEFFICIENT
            Y = pdist(distance_matrix_man)
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
            linkage_matrix = shc.linkage(distance_matrix_man, method =  cluster[i])
            shc.set_link_color_palette(['#d81159' ,  '#8f2d56' ,   '#218380' ,   '#fbb13c' ,   '#73d2de'])

            # COPHENETIC COEFFICIENT
            Y = pdist(distance_matrix_man)
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
            linkage_matrix = shc.linkage(distance_matrix_man, method =  cluster[i])
            shc.set_link_color_palette(['#d81159' ,  '#8f2d56' ,   '#218380' ,   '#fbb13c' ,   '#73d2de'])

            # COPHENETIC COEFFICIENT
            Y = pdist(distance_matrix_man)
            Z = linkage_matrix
            c, coph_dists = cophenet(Z, Y)
            print(f"Cophenetic correlation coefficient for {cluster[i]}: {c}")
        

            plt.subplot(2,3,i+1)
            plt.title("Hierarquy Cluster:  " + " " + cluster[i], fontsize = 10, color = 'black')
            dendrogram = shc.dendrogram(linkage_matrix, orientation='left', labels =labels)
               

    plt.suptitle("Manhattan Distance")
    plt.subplots_adjust(left=0.1, right=0.9, bottom=0.1, top=0.9, wspace=0.5, hspace=0.5)
    plt.subplots.figsize = (20,40)
    plt.show()

# ------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------

# NORMALIZED MANHATTAN DISTANCE ------------------------------------------------------------------------------------------------------------

def word_based_method_man_norm(seq1,seq2,k):
    """
    Given two sequence-representing vectors, returns the Manhattan normalized distance between them.
    """
    
    w1 = f.generate_kmers(seq1, k)
    w2 = f.generate_kmers(seq2, k)
    
    w_union = f.union_vectors(w1, w2)

    count1= f.count_items_norm(w_union, w1)
    count2= f.count_items_norm(w_union, w2)

    distance = manhattan(count1,count2)
    
    return distance

def best_k_man_norm(seq1,seq2):
    
    """
    Given two sequence-representing vectors, returns the k that minimizes the Manhattan normalized distance between them (and the Bray Curtis distance)
    """   
    
    kmax=15
    dist = word_based_method_man_norm(seq1,seq2,1)
    k=1
    
    for i in range(2,kmax + 1):
        if dist >word_based_method_man_norm(seq1,seq2,i):
            k = i
            dist = word_based_method_man_norm(seq1,seq2,i)

    return [k,dist]

def all_k_man_norm(seq1,seq2):
    kmax=min(len(seq1),len(seq2))
    dist = []
    k =list(range(1,kmax+1))
    
    for i in range(1,kmax+1):
        dist += [word_based_method_man_norm(seq1,seq2,i)]
    
# NORMALIZED MANHATTAN DISTANCE FOR MORE THAN 2 SEQUENCES

def matrix_k_man_norm(vector):
    
    """
    Given a vector with n sequences,creates a matrix where M(a,b) is the best k for the normalized Manhattan distance
    """
    
    n =len(vector)
    M = np.zeros((n,n))
    for i in range(n):
        for j in range(i,n):
            M[i,j] = best_k_man_norm(vector[i],vector[j])[0]
            
    median = int(np.median(M))
            
    return [median,M]

def matrix_dist_man_norm(vector):
    k = matrix_k_man_norm(vector)[0]
    n =len(vector)
    M = np.zeros((n,n))
    
    for i in range(n):
        for j in range(n):
            M[i,j] = word_based_method_man_norm(vector[i],vector[j],k)
            
    return M


def plot_man_norm(genes,labels,cluster):
    n_plots = len(cluster)
    distance_matrix_man = matrix_dist_man_norm(genes)

    if n_plots <=3:
        for i in range(n_plots):
            # calculate the linkage matrix
            linkage_matrix = shc.linkage(distance_matrix_man, method =  cluster[i])
            shc.set_link_color_palette(['#d81159' ,  '#8f2d56' ,   '#218380' ,   '#fbb13c' ,   '#73d2de'])

            # COPHENETIC COEFFICIENT
            Y = pdist(distance_matrix_man)
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
            linkage_matrix = shc.linkage(distance_matrix_man, method =  cluster[i])
            shc.set_link_color_palette(['#d81159' ,  '#8f2d56' ,   '#218380' ,   '#fbb13c' ,   '#73d2de'])

            # COPHENETIC COEFFICIENT
            Y = pdist(distance_matrix_man)
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
            linkage_matrix = shc.linkage(distance_matrix_man, method =  cluster[i])
            shc.set_link_color_palette(['#d81159' ,  '#8f2d56' ,   '#218380' ,   '#fbb13c' ,   '#73d2de'])
        
            # COPHENETIC COEFFICIENT
            Y = pdist(distance_matrix_man)
            Z = linkage_matrix
            c, coph_dists = cophenet(Z, Y)
            print(f"Cophenetic correlation coefficient for {cluster[i]}: {c}")

            plt.subplot(2,3,i+1)
            plt.title("Hierarquy Cluster:  " + " " + cluster[i], fontsize = 10, color = 'black')
            dendrogram = shc.dendrogram(linkage_matrix, orientation='left', labels =labels)
               

    plt.suptitle("Normalized Manhattan Distance")
    plt.subplots_adjust(left=0.1, right=0.9, bottom=0.1, top=0.9, wspace=0.5, hspace=0.5)
    plt.subplots.figsize = (20,40)
    plt.show()

#method = ["complete" , "single"]
#Labels =  ['epsilon', 'G_gamma', 'A_gamma', 'delta', 'beta']
#Sequences = ['MVHFTAEEKAAVTSLWSKMNVEEAGGEALGRLLVVYPWTQRFFDSFGNLSSPSAILGNPKVKAHGKKVLTSFGDAIKNMDNLKPAFAKLSELHCDKLHVDPENFKLLGNVMVIILATHFGKEFTPEVQAAWQKLVSAVAIALAHKYH', 'MGHFTEEDKATITSLWGKVNVEDAGGETLGRLLVVYPWTQRFFDSFGNLSSASAIMGNPKVKAHGKKVLTSLGDAIKHLDDLKGTFAQLSELHCDKLHVDPENFKLLGNVLVTVLAIHFGKEFTPEVQASWQKMVTGVASALSSRYH', 'MGHFTEEDKATITSLWGKVNVEDAGGETLGRLLVVYPWTQRFFDSFGNLSSASAIMGNPKVKAHGKKVLTSLGDAIKHLDDLKGTFAQLSELHCDKLHVDPENFKLLGNVLVTVLAIHFGKEFTPEVQASWQKMVTAVASALSSRYH', 'MVHLTPEEKTAVNALWGKVNVDAVGGEALGRLLVVYPWTQRFFESFGDLSSPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFSQLSELHCDKLHVDPENFRLLGNVLVCVLARNFGKEFTPQMQAAYQKVVAGVANALAHKYH', 'MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH']
#plot_man(Sequences, Labels, method)

