# This is a sample Python script.
from Data_PCA import exp_pca
from Data_GAE import *
from Data_CoReg import *


if __name__ == '__main__':
    # Gene expression data processing
    exp_pca()

    # The gene weights are calculated based on the PPI network
    Running_ppi_gae()

    #The gene weights are calculated based on the PPI network

    ##The regulatory relationships emerging from PPI networks were screened
    reg_save()

    ##The gene weights are calculated according to PageRank
    pageRank_cal()

    ##Co-regulatory values are calculated based on regulatory relationships
    coRegCal()