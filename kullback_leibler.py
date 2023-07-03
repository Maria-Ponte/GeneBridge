import funcoes_gerais as fg
from scipy.special import rel_entr
import numpy as np
import math
import scipy.cluster.hierarchy as shc
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import single, cophenet
from scipy.spatial.distance import pdist, squareform


x = "ATGTGTG"
y = "CATGTG"

def kullback_leibler(seq1, seq2,k): 
    w1 = sorted(list(set(fg.generate_kmers(seq1,k))))
    w2 = sorted(list(set(fg.generate_kmers(seq2,k))))
    w12 = sorted(list(set(w1+w2)))
    
    p1 = fg.count_items_norm(w12,seq1)
    p2 = fg.count_items_norm(w12,seq2)

    kl = 0
    for i in range(len(w12)):
        if p1[i]*p2[i] != 0:
            kl += p1[i]*math.log(p1[i]/p2[i])
    return sum(rel_entr(p1,p2))


def distance_matrix(sequences,k):

    num_sequences = len(sequences)
    distance_matrix = np.zeros((num_sequences, num_sequences))

    for i in range(num_sequences):
        for j in range(num_sequences):
            distance_matrix[i][j] = kullback_leibler(sequences[i], sequences[j],k)

    return distance_matrix

def infinite_values(matrix):
   
    if np.all(np.isfinite(matrix)):
        matrix = matrix
    
    else: 
        max_value = np.max(matrix[np.isfinite(matrix)])
        matrix[np.isinf(matrix)] = max_value*2

    return  matrix

def plot_kl(genes,labels, cluster):
    
    n_plots = len(cluster)
    dm = infinite_values(distance_matrix(genes , 1))           #ESTOU A ASSUMIR K=1 MAS SE CALHAR PODEMOS MUDAR ISTO
    print(dm)

    if n_plots <=3:
        for i in range(n_plots):
            # calculate the linkage matrix
            linkage_matrix = shc.linkage(dm, method =  cluster[i])
            shc.set_link_color_palette(['#d81159' ,  '#8f2d56' ,   '#218380' ,   '#fbb13c' ,   '#73d2de'])

            # COPHENETIC COEFFICIENT
            Y = pdist(dm)
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
            linkage_matrix = shc.linkage(dm, method =  cluster[i])
            shc.set_link_color_palette(['#d81159' ,  '#8f2d56' ,   '#218380' ,   '#fbb13c' ,   '#73d2de'])

            # COPHENETIC COEFFICIENT
            Y = pdist(dm)
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
            linkage_matrix = shc.linkage(dm, method =  cluster[i])
            shc.set_link_color_palette(['#d81159' ,  '#8f2d56' ,   '#218380' ,   '#fbb13c' ,   '#73d2de'])
            
            # COPHENETIC COEFFICIENT
            Y = pdist(dm)
            Z = linkage_matrix
            c, coph_dists = cophenet(Z, Y)
            print(f"Cophenetic correlation coefficient for {cluster[i]}: {c}")

            plt.subplot(2,3,i+1)
            plt.title("Hierarquy Cluster:  " + " " + cluster[i], fontsize = 10, color = 'black')
            dendrogram = shc.dendrogram(linkage_matrix, orientation='left', labels =labels)
               

    plt.suptitle("Kullback-Leibler Divergence")
    plt.subplots_adjust(left=0.1, right=0.9, bottom=0.1, top=0.9, wspace=0.5, hspace=0.5)
    plt.subplots.figsize = (20,40)
    plt.show()
    

#method = ["complete", "single", "average"]
#Labels =  ['epsilon', 'G_gamma', 'A_gamma', 'delta', 'beta']
#Sequences = ['MVHFTAEEKAAVTSLWSKMNVEEAGGEALGRLLVVYPWTQRFFDSFGNLSSPSAILGNPKVKAHGKKVLTSFGDAIKNMDNLKPAFAKLSELHCDKLHVDPENFKLLGNVMVIILATHFGKEFTPEVQAAWQKLVSAVAIALAHKYH', 'MGHFTEEDKATITSLWGKVNVEDAGGETLGRLLVVYPWTQRFFDSFGNLSSASAIMGNPKVKAHGKKVLTSLGDAIKHLDDLKGTFAQLSELHCDKLHVDPENFKLLGNVLVTVLAIHFGKEFTPEVQASWQKMVTGVASALSSRYH', 'MGHFTEEDKATITSLWGKVNVEDAGGETLGRLLVVYPWTQRFFDSFGNLSSASAIMGNPKVKAHGKKVLTSLGDAIKHLDDLKGTFAQLSELHCDKLHVDPENFKLLGNVLVTVLAIHFGKEFTPEVQASWQKMVTAVASALSSRYH', 'MVHLTPEEKTAVNALWGKVNVDAVGGEALGRLLVVYPWTQRFFESFGDLSSPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFSQLSELHCDKLHVDPENFRLLGNVLVCVLARNFGKEFTPQMQAAYQKVVAGVANALAHKYH', 'MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH']
