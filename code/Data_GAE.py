import numpy as np
import scipy.sparse as sp
import torch
import math
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
import random
import time
from torch.optim import Adam
from sklearn import preprocessing


def encode_onehot(labels):
    classes = set(labels)
    class_dict = {c: np.identity(len(classes))[i, :] for i, c in enumerate(classes)}

    label_onehot = np.array(list(map(class_dict.get, labels)),
                            dtype=np.int32)

    return label_onehot


def normalize(mx):
    rowsum = np.array(mx.sum(1))
    r_inv = np.power(rowsum, -1).flatten()
    r_inv[np.isinf(r_inv)] = 0
    r_mat_inv = sp.diags(r_inv)
    mx = r_mat_inv.dot(mx)

    return mx


def normalize_adj(adjacency):
    degree = np.array(adjacency.sum(1))
    d_hat = sp.diags(np.power(degree, -0.5).flatten())
    adj_norm = d_hat.dot(adjacency).dot(d_hat).tocoo()
    return adj_norm


def sparse_to_tuple(sparse_mx):
    if not sp.isspmatrix_coo(sparse_mx):  # 判断是否为csr_matrix类型
        sparse_mx = sparse_mx.tocoo()  # 实现类型转换
    coords = np.vstack((sparse_mx.row, sparse_mx.col)).transpose()
    # np.vstack按垂直方向（行顺序）堆叠数组构成一个新的数组，堆叠的数组需要具有相同的维度，transpose()作用是转置
    values = sparse_mx.data
    shape = sparse_mx.shape
    return coords, values, shape


def l1_penalty(w):
    return (torch.abs(w)).sum()


def l2_penalty(w):
    return (w ** 2).sum() / 2


def sparse_mx_to_torch_sparse_tensor(sparse_mx):
    sparse_mx = sparse_mx.tocoo().astype(np.float32)
    indices = torch.from_numpy(
        np.vstack((sparse_mx.row, sparse_mx.col)).astype(np.int64))
    values = torch.from_numpy(sparse_mx.data)
    shape = torch.Size(sparse_mx.shape)

    return torch.sparse.FloatTensor(indices, values, shape)


def load_data(data1, data2):
    print("Loading data...")

    idx_features_labels = np.genfromtxt(data1, delimiter=",", dtype=np.dtype(str))
    features = sp.csr_matrix(idx_features_labels[:, 1:-1], dtype=np.float32)
    labels = encode_onehot(idx_features_labels[:, -1])

    # build graph
    idx = np.array(idx_features_labels[:, 0], dtype=np.int32)
    idx_map = {j: i for i, j in enumerate(idx)}
    edges_unordered = np.genfromtxt(data2, delimiter="\t", dtype=np.int32)
    edges = np.array(list(map(idx_map.get, edges_unordered.flatten())), dtype=np.int32,
                     ).reshape(edges_unordered.shape)
    adj = sp.coo_matrix((np.ones(edges.shape[0]), (edges[:, 0], edges[:, 1])),
                        shape=(labels.shape[0], labels.shape[0]),
                        dtype=np.float32)

    adj = adj + adj.T.multiply(adj.T > adj) - adj.multiply(adj.T > adj)


    adj = normalize_adj(adj + sp.eye(adj.shape[0]))

    features = torch.FloatTensor(np.array(features.todense()))
    labels = torch.LongTensor(np.where(labels)[1])
    adj = sparse_mx_to_torch_sparse_tensor(adj)

    print("Data loaded successfully...")

    return adj, features, labels


adj, features, labels = load_data(data1="data/gae_RNA_MaxMin.csv",
                                  data2="data/gae_edges.txt")


class GraphConvSparse(nn.Module):
    def __init__(self, input_dim, output_dim, adj, activation=F.leaky_relu):
        super(GraphConvSparse, self).__init__()
        self.weight = nn.Parameter(torch.FloatTensor(input_dim, output_dim))
        self.adj = adj
        self.reset_parameters()
        self.activation = activation

    def reset_parameters(self):  # 重置参数
        nn.init.kaiming_uniform_(self.weight)

    def forward(self, inputs):
        x = inputs
        x = torch.mm(x, self.weight)  # x与系数矩阵相乘
        x = torch.mm(self.adj, x)  # 邻接矩阵与x相乘
        outputs = self.activation(x)
        return outputs


def dot_product_decode(Z):
    A_pred = torch.sigmoid(torch.matmul(Z, Z.t()))
    # A_pred = torch.matmul(Z, Z.t())
    return A_pred


class GAE(nn.Module):
    def __init__(self, adj):
        super(GAE, self).__init__()
        self.base_gcn = GraphConvSparse(input_dim, hidden1_dim, adj)
        self.gcn_mean = GraphConvSparse(hidden1_dim, hidden2_dim, adj, activation=lambda x: x)

        # zsp
        self.decoder_base_gcn = GraphConvSparse(hidden2_dim, hidden1_dim, adj)
        self.decoder_gcn_mean = GraphConvSparse(hidden1_dim, input_dim, adj, activation=lambda x: x)

    def encode(self, X):
        hidden = self.base_gcn(X)
        z = self.mean = self.gcn_mean(hidden)
        return z

    # zsp
    def decode(self, z):
        hidden = self.decoder_base_gcn(z)
        z = self.decoder_gcn_mean(hidden)
        return z

    def forward(self, X):
        Z = self.encode(X)
        out_file = Z
        result = self.decode(Z)
        A_pred = dot_product_decode(Z)
        return result, A_pred, out_file

def Running_ppi_gae():

    adj, features, labels = load_data(data1="data/gae_RNA_MaxMin.csv",
                                      data2="data/gae_edges.txt")

    input_dim = 789
    hidden1_dim = 200
    hidden2_dim = 50
    use_feature = True
    num_epoch = 200
    learning_rate = 0.0001

    model = GAE(adj)
    optimizer = Adam(model.parameters(), lr=learning_rate)

    torch.save(model.state_dict(), 'data/gae.pth')

    for epoch in range(num_epoch):
        t = time.time()

        feature_pre, A_pred, output = model(features)
        optimizer.zero_grad()

        loss_struct = F.binary_cross_entropy(A_pred.view(-1), adj.to_dense().view(-1))

        loss_feature = F.mse_loss(feature_pre, features)

        loss = loss_feature + loss_struct

        loss.backward()
        optimizer.step()

        print("Epoch:", '%04d' % (epoch + 1),
              "train_loss=", "{:.5f}".format(loss.item()),
              "features_loss=", "{:.5f}".format(loss_feature),
              "struct_loss=", "{:.5f}".format(loss_struct),
              "time=", "{:.5f}".format(time.time() - t))

        if ((epoch + 1) % 5 == 0):
            print(output.detach().numpy().shape)
            np.savetxt('data/ppi_gae_%s_%s.csv' % ((epoch + 1), loss.data), output.detach().numpy(),
                       delimiter=',')
