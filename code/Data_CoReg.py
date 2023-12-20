import networkx as nx
import numpy as np


# reg dict
def reg_save():
    gene_set = set()
    trrust = np.loadtxt("data/TRRUST.txt", delimiter='\t', dtype=str)
    gene_id = np.loadtxt('data/largest_component.csv', delimiter='\t', dtype=str)

    for i in trrust:
        for j in i:
            if j in gene_id:
                gene_set.add(j)

    print(len(gene_set))
    np.savetxt("data/coReg/reg_name.csv", np.array(list(gene_set)), delimiter=',', fmt='%s')

    reg_dict = dict()
    for gene in gene_set:
        print(gene)
        a_list = []
        for i in range(trrust.shape[0]):
            if trrust[i][1] == gene and trrust[i][0] in gene_id:
                a_list.append(trrust[i][0])
            reg_dict[gene] = a_list
        # if len(a_list) == 0:
        #     print(gene)
    # print(reg_dict)

    print(len(reg_dict))
    np.save('data/coReg/reg_dict.npy', reg_dict)


# reg_save()


def pageRank_cal():
    edges_file = np.loadtxt("data/largest_edges.csv", delimiter=',', dtype=str)
    print(str(edges_file.shape))

    G = nx.DiGraph()
    for i in range(0, edges_file.shape[0]):
        G.add_edge(edges_file[i][0], edges_file[i][1])
        G.add_edge(edges_file[i][1], edges_file[i][0])

    pageRank_list = nx.pagerank(G, alpha=1)
    print(len(pageRank_list))
    print(type(pageRank_list))

    pageRank_order = sorted(pageRank_list.items(), key=lambda x: x[1], reverse=True)
    np.savetxt("data/coReg/nodePageRank.csv", np.array(pageRank_order), delimiter=',', fmt='%s')

    reg_name = np.loadtxt('data/coReg/reg_name.csv', delimiter=',', dtype=str)
    p6 = []

    for i in range(len(pageRank_order)):
        print(i)
        if pageRank_order[i][0] in reg_name:
            p6.append(pageRank_order[i][0])

    print(len(p6))

    weight = []
    for j in range(0, len(p6)):
        weight.append([p6[j], 1 / (j + 1) ** (1 / 2)])
    print(len(weight))
    # print(weight)
    np.savetxt('data/coReg/node_Rank.csv', np.array(weight), delimiter=',', fmt='%s')


# pageRank_cal()


def coRegCal():
    b_set = np.loadtxt('data/coReg/reg_name.csv', delimiter='\t', dtype=str)
    # b_set = np.loadtxt('data/reg_name.txt',delimiter='\t',dtype=str)
    print(b_set.shape)

    load_dict = np.load('data/coReg/reg_dict.npy', allow_pickle=True).item()

    node_pagerank = np.loadtxt('data/coReg/node_Rank.csv', delimiter=',', dtype=str)
    print(node_pagerank.shape)
    # b_set = []
    # for i in reg_name:
    #     if i in node_pagerank[:, 0]:
    #         b_set.append(i)

    print("b_set的长度：" + str(len(b_set)))

    coReg_list = []
    for i in range(b_set.shape[0]):
        print(i)
        a = load_dict[b_set[i]]
        a_rank = float(node_pagerank[(np.where(node_pagerank[:, 0] == b_set[i])[0]), 1])

        a_weight = []
        for a_item in a:
            for n in range(node_pagerank.shape[0]):
                if node_pagerank[n][0] == a_item:
                    # a_node.append(node_pagerank[n][0])
                    a_weight.append(float(node_pagerank[n][1]))
        for j in range(i, len(b_set)):
            b = load_dict[b_set[j]]
            b_rank = float(node_pagerank[(np.where(node_pagerank[:, 0] == b_set[j])[0]), 1])
            b_weight = []
            c = list(set(a).intersection(set(b)))
            # c_node = []
            c_weight = []
            if b_set[i] != b_set[j] and len(c) != 0:
                for b_item in b:
                    for n in range(node_pagerank.shape[0]):
                        if node_pagerank[n][0] == b_item:
                            # c_node.append(node_pagerank[n][0])
                            b_weight.append(float(node_pagerank[n][1]))
                for m in c:
                    for n in range(node_pagerank.shape[0]):
                        if node_pagerank[n][0] == m:
                            # c_node.append(node_pagerank[n][0])
                            c_weight.append(float(node_pagerank[n][1]))

                if len(c_weight) != 0:
                    coReg_list.append([b_set[i], b_set[j],2 * sum(c_weight) / (sum(a_weight) + sum(b_weight)) * (a_rank * b_rank) ** (1 / 2)])
    # print(coReg_list)
    print(len(coReg_list))
    np.savetxt('data/node_coReg.txt', np.array(coReg_list), delimiter='\t', fmt='%s')

# coRegCal()
