import numpy as np
import math



def load_modul_list():
    filename = 'output/com_p10.txt'
    modul_list = []
    score = np.loadtxt('output/score_p10.txt', delimiter='\t', dtype=str)
    with open(filename) as f:
        lines = f.readlines()
        for line in lines:
            line = line.split()
            modul_list.append(line)
    print(len(modul_list))
    print(modul_list)

    score_1 = []  # 直接均值
    score_2 = []  # 考虑size
    for i in range(len(modul_list)):
        list_size = len(modul_list[i])
        # print(list_size)
        list_edges = len(modul_list[i]) * (len(modul_list[i]) - 1) / 2
        list_score = float(score[i])

        current_score1 = list_score / list_edges
        x = math.log(1 + 1 / (1 - current_score1 ** 2))
        s = (list_size / 2) ** (1 / 2)
        z = 4 * math.exp(s) / (1 + 4 * math.exp(s))
        # current_score2 = (list_size/2)**2*0.5*x
        current_score2 = z * x

        score_1.append(current_score1)
        score_2.append(current_score2)
        # score_2.append((list_num ** (1 / 4)) * list_score / (2 * list_edges))
        modul_list[i].append(str(current_score1))

        # modul_list[i].append(str(x))
        #
        # modul_list[i].append(str(z))

        modul_list[i].append(str(current_score2))
        # modul_list[i].append((list_num ** (1 /4)) * list_score / (2 * list_edges))
        # print(list_num)

    score1 = np.array(score_1)
    s1 = (-score1).argsort()
    score2 = np.array(score_2)
    s2 = (-score2).argsort()
    modul_list_score1 = []
    modul_list_score2 = []
    for i in s1:
        for j in range(len(modul_list)):
            if j == i:
                modul_list_num = modul_list[j]
                modul_list_score1.append(modul_list_num)
    for i in s2:
        for j in range(len(modul_list)):
            if j == i:
                modul_list_num = modul_list[j]
                modul_list_score2.append(modul_list_num)

    print(score_1)
    print(score_2)
    print(len(score_1))
    print(len(score_2))

    save_list('output/p10_4_DMMO.txt', modul_list_score2)


load_modul_list()

