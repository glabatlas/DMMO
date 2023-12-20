import numpy as np
from sklearn.decomposition import PCA

def exp_pca():
    # load the gene expression data
    exp_data = np.loadtxt('data/RNA_exp.csv', dtype=str, delimiter=',')
    print("Data dimension： " + str(exp_data.shape))

    # load the gene name
    genes_id = exp_data[:, 0]
    genes_id = genes_id.reshape(-1, 1)
    print("Number of genes： " + str(genes_id.shape))

    # load data
    exp_data = np.delete(exp_data, 0, axis=1)
    exp_data = exp_data.astype(float)
    print("data dimension“ " + str(exp_data.shape))


    # model
    pca = PCA(n_components = 10)

    exp_pca = pca.fit_transform(exp_data)
    print(exp_pca.shape)
    print(pca.explained_variance_ratio_)
    print("Variance contribution rate:  " + str(sum(pca.explained_variance_ratio_)))

    print(exp_pca.shape)
    exp_pca_data = np.concatenate((genes_id, exp_pca), axis=1)  # dimension-reduced part of the data is merged with the gene name
    print(exp_pca_data.shape)
    np.savetxt('output/exp_pca.csv', exp_pca_data, delimiter=',', fmt='%s')


# exp_pca()
