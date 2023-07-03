
import numpy as np
import funcoes_gerais as f
import scipy.cluster.hierarchy as shc
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import fcluster
from sklearn.metrics import silhouette_score
from scipy.cluster.hierarchy import single, cophenet
from scipy.spatial.distance import pdist, squareform




def euclidean_distance(vector1, vector2):
    """
    Given two vectors, return the Euclidean distance between them.
    """
    vector1 = np.array(vector1)
    vector2 = np.array(vector2)
    distance = np.linalg.norm(vector1 - vector2)
    return distance


# EUCLIDIAN DISTANCE -----------------------------------------------------------------------------------------------------------------------

def word_based_method_eucl(seq1,seq2,k):
    """
    Given two sequence-representing vectors, returns the Euclidean distance between them.
    """
    
    w1 = f.generate_kmers(seq1, k)
    w2 = f.generate_kmers(seq2, k)
    
    w_union = f.union_vectors(w1, w2)
    
    count1= f.count_items(w_union, w1)
    count2= f.count_items(w_union, w2)

    distance = euclidean_distance(count1,count2)
    
    return distance

def best_k_eucl(seq1,seq2):
    
    """
    Given two sequence-representing vectors, returns the k that minimizes the Euclidean distance between them (and the Euclidean distance)
    """   
    
    kmax=15
    dist = word_based_method_eucl(seq1,seq2,1)
    k=1
    
    for i in range(2,kmax + 1):
        if dist >word_based_method_eucl(seq1,seq2,i):
            k = i
            dist = word_based_method_eucl(seq1,seq2,i)

    return k,dist

def all_k_eucl(seq1,seq2):
    kmax=min(len(seq1),len(seq2))
    dist = []
    k =list(range(1,kmax+1))
    
    for i in range(1,kmax+1):
        dist += [word_based_method_eucl(seq1,seq2,i)]
    
    return k,dist

# EUCLIDIAN DISTANCE FOR MORE THAN 2 SEQUENCES

def matrix_k_eucl(vector):
    
    """
    Given a vector with n sequences,creates a matrix where M(a,b) is the best k for the normalized eucladian distance
    """
    
    n =len(vector)
    M = np.zeros((n,n))
    for i in range(n):
        for j in range(i,n):
            M[i,j] = best_k_eucl(vector[i],vector[j])[0]
            
    median = int(np.median(M))
            
    return [median,M]

def matrix_dist_eucl(vector):
    k = matrix_k_eucl(vector)[0]
    n =len(vector)
    M = np.zeros((n,n))
    
    for i in range(n):
        for j in range(n):
            M[i,j] = word_based_method_eucl(vector[i],vector[j],k)
            
    return M

def plot_eucl(genes,labels, cluster):
    n_plots = len(cluster)
    
    distance_matrix_eucl = matrix_dist_eucl(genes)
    
    if n_plots <=3:
        for i in range(n_plots):
            # calculate the linkage matrix
            linkage_matrix = shc.linkage(distance_matrix_eucl, method =  cluster[i])
            shc.set_link_color_palette(['#d81159' ,  '#8f2d56' ,   '#218380' ,   '#fbb13c' ,   '#73d2de'])

            # COPHENETIC COEFFICIENT
            Y = pdist(distance_matrix_eucl)
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
            linkage_matrix = shc.linkage(distance_matrix_eucl, method =  cluster[i])
            shc.set_link_color_palette(['#d81159' ,  '#8f2d56' ,   '#218380' ,   '#fbb13c' ,   '#73d2de'])

            # COPHENETIC COEFFICIENT
            Y = pdist(distance_matrix_eucl)
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
            linkage_matrix = shc.linkage(distance_matrix_eucl, method =  cluster[i])
            shc.set_link_color_palette(['#d81159' ,  '#8f2d56' ,   '#218380' ,   '#fbb13c' ,   '#73d2de'])
            
            # COPHENETIC COEFFICIENT
            Y = pdist(distance_matrix_eucl)
            Z = linkage_matrix
            c, coph_dists = cophenet(Z, Y)
            print(f"Cophenetic correlation coefficient for {cluster[i]}: {c}")
            
            plt.subplot(2,3,i+1)
            plt.title("Hierarquy Cluster:  " + " " + cluster[i], fontsize = 10, color = 'black')
            dendrogram = shc.dendrogram(linkage_matrix, orientation='left', labels =labels)

    

    plt.suptitle("Euclidean Distance")
    plt.subplots_adjust(left=0.1, right=0.9, bottom=0.1, top=0.9, wspace=0.5, hspace=0.5)
    plt.subplots.figsize = (20,40)
    plt.show()

# ------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------

# NORMALIZED EUCLIDIAN DISTANCE ------------------------------------------------------------------------------------------------------------

def word_based_method_eucl_norm(seq1,seq2,k):
    """
    Given two sequence-representing vectors, returns the Euclidean distance between them.
    """
    
    w1 = f.generate_kmers(seq1, k)
    w2 = f.generate_kmers(seq2, k)
    
    w_union = f.union_vectors(w1, w2)

    count1= f.count_items_norm(w_union, w1)
    count2= f.count_items_norm(w_union, w2)

    distance = euclidean_distance(count1,count2)
    
    return distance

def best_k_eucl_norm(seq1,seq2):
    
    """
    Given two sequence-representing vectors, returns the k that minimizes the Euclidean distance between them (and the Euclidean distance)
    """   
    
    kmax=15
    dist = word_based_method_eucl_norm(seq1,seq2,1)
    k=1
    
    for i in range(2,kmax + 1):
        if dist >word_based_method_eucl_norm(seq1,seq2,i):
            k = i
            dist = word_based_method_eucl_norm(seq1,seq2,i)

    return [k,dist]

def all_k_eucl_norm(seq1,seq2):
    kmax=min(len(seq1),len(seq2))
    dist = []
    k =list(range(1,kmax+1))
    
    for i in range(1,kmax+1):
        dist += [word_based_method_eucl_norm(seq1,seq2,i)]
    
    return k,dist

# NORMALIZED EUCLIDIAN DISTANCE FOR MORE THAN 2 SEQUENCES

def matrix_k_eucl_norm(vector):
    
    """
    Given a vector with n sequences,creates a matrix where M(a,b) is the best k for the normalized eucladian distance
    """
    
    n =len(vector)
    M = np.zeros((n,n))
    for i in range(n):
        for j in range(i,n):
            M[i,j] = best_k_eucl_norm(vector[i],vector[j])[0]
            
    median = int(np.median(M))
            
    return [median,M]

def matrix_dist_eucl_norm(vector):
    k = matrix_k_eucl_norm(vector)[0]
    n =len(vector)
    M = np.zeros((n,n))
    
    for i in range(n):
        for j in range(n):
            M[i,j] = word_based_method_eucl_norm(vector[i],vector[j],k)
            
    return M

def plot_eucl_norm(genes,labels, cluster):
    n_plots = len(cluster)
    
    distance_matrix_eucl = matrix_dist_eucl_norm(genes)

    
    if n_plots <=3:
        for i in range(n_plots):
            # calculate the linkage matrix
            linkage_matrix = shc.linkage(distance_matrix_eucl, method =  cluster[i])
            shc.set_link_color_palette(['#d81159' ,  '#8f2d56' ,   '#218380' ,   '#fbb13c' ,   '#73d2de'])
           
            # COPHENETIC COEFFICIENT
            Y = pdist(distance_matrix_eucl)
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
            linkage_matrix = shc.linkage(distance_matrix_eucl, method =  cluster[i])
            shc.set_link_color_palette(['#d81159' ,  '#8f2d56' ,   '#218380' ,   '#fbb13c' ,   '#73d2de'])

            # COPHENETIC COEFFICIENT
            Y = pdist(distance_matrix_eucl)
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
            linkage_matrix = shc.linkage(distance_matrix_eucl, method =  cluster[i])
            shc.set_link_color_palette(['#d81159' ,  '#8f2d56' ,   '#218380' ,   '#fbb13c' ,   '#73d2de'])

            # COPHENETIC COEFFICIENT
            Y = pdist(distance_matrix_eucl)
            Z = linkage_matrix
            c, coph_dists = cophenet(Z, Y)
            print(f"Cophenetic correlation coefficient for {cluster[i]}: {c}")
        

            plt.subplot(2,3,i+1)
            plt.title("Hierarquy Cluster:  " + " " + cluster[i], fontsize = 10, color = 'black')
            dendrogram = shc.dendrogram(linkage_matrix, orientation='left', labels =labels)
               

    plt.suptitle("Normalized Euclidean Distance")
    plt.subplots_adjust(left=0.1, right=0.9, bottom=0.1, top=0.9, wspace=0.5, hspace=0.5)
    plt.subplots.figsize = (20,40)
    plt.show()

    

#method = ["complete"]
#Labels =  ['epsilon', 'G_gamma', 'A_gamma', 'delta', 'beta']
#Sequences = ['MVHFTAEEKAAVTSLWSKMNVEEAGGEALGRLLVVYPWTQRFFDSFGNLSSPSAILGNPKVKAHGKKVLTSFGDAIKNMDNLKPAFAKLSELHCDKLHVDPENFKLLGNVMVIILATHFGKEFTPEVQAAWQKLVSAVAIALAHKYH', 'MGHFTEEDKATITSLWGKVNVEDAGGETLGRLLVVYPWTQRFFDSFGNLSSASAIMGNPKVKAHGKKVLTSLGDAIKHLDDLKGTFAQLSELHCDKLHVDPENFKLLGNVLVTVLAIHFGKEFTPEVQASWQKMVTGVASALSSRYH', 'MGHFTEEDKATITSLWGKVNVEDAGGETLGRLLVVYPWTQRFFDSFGNLSSASAIMGNPKVKAHGKKVLTSLGDAIKHLDDLKGTFAQLSELHCDKLHVDPENFKLLGNVLVTVLAIHFGKEFTPEVQASWQKMVTAVASALSSRYH', 'MVHLTPEEKTAVNALWGKVNVDAVGGEALGRLLVVYPWTQRFFESFGDLSSPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFSQLSELHCDKLHVDPENFRLLGNVLVCVLARNFGKEFTPQMQAAYQKVVAGVANALAHKYH', 'MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH']

