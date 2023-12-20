library(parallel)
library(plyr)

normalize<-function(scores)
{
  
  # Z-score
  mean_scores <- mean(scores)
  var_scores <- var(scores)
  
  if (var_scores == 0) {
    var_scores <- 1
    mean_scores <- 0
  }
  
  zScores <- (scores-mean_scores)/sqrt(var_scores)
  
  # Scale to  (0,1)
  if (min(zScores) != max(zScores)) {
    min_zScores <- min(zScores)
    zScores <- zScores - min_zScores
    max_zScores <- max(zScores)
    if (max_zScores != 0) {zScores <- zScores / max_zScores}
  }
  
  return (zScores)
}

mut_function <- function(gene1){
  out_result <- vector(length = noGene)
  print(paste("current gene mut:", gene1))
  for (gene2 in (gene1+1):noGene){
    cur <- mut[,c(gene1,gene2)]
    mutScore <- 2 * sum(rowSums(cur) != 0) - sum(cur)
    out_result[gene2] <- mutScore 
  }
  return(out_result)
}

k_3 <- function(rank_num){
  gene1 <- k_2_res[rank_num,1]
  gene2 <- k_2_res[rank_num,2]
  # print(paste("current gene1 :", gene1))
  current_score <- matrix(data = 0,ncol = 4,nrow = 2000)
  print(paste("  current time  ",Sys.time(),"current K3 rank  ", rank_num))
  # print(paste("current gene :", gene1,"  -  ",gene2," -  ",gene3))
  for (gene3 in 1:noGene) {
    if(gene3!=gene1&&gene3!=gene2){
      score <- sum(pairsZScore[c(gene1,gene2,gene3),c(gene1,gene2,gene3)]) 
      if(score > min(current_score[,4])){
        current_score[which.min(current_score[,4]),1:3] = sort(c(gene1,gene2,gene3))
        current_score[which.min(current_score[,4]),4] = score
      }
    }
  }
  return(current_score)
}

k_4 <- function(rank_num){
  gene1 <- k_3_res[rank_num,1]
  gene2 <- k_3_res[rank_num,2]
  gene3 <- k_3_res[rank_num,3]
  # print(paste("current gene1 :", gene1))
  current_score <- matrix(data = 0,ncol = 5,nrow = 2000)
  print(paste("  current time:",Sys.time(),"current k4 rank:", rank_num))
  # print(paste("current gene :", gene1,"  -  ",gene2," -  ",gene3))
  for (gene4 in 1:noGene) {
    if(gene4!=gene1&&gene4!=gene2&&gene4!=gene3){
      score <- sum(pairsZScore[c(gene1,gene2,gene3,gene4),c(gene1,gene2,gene3,gene4)])
      if(score > min(current_score[,5])){
        current_score[which.min(current_score[,5]),1:4] = sort(c(gene1,gene2,gene3,gene4))
        current_score[which.min(current_score[,5]),5] = score
      }
    }
  }
  return(current_score)
}

k_5 <- function(rank_num){
  gene1 <- k_4_res[rank_num,1]
  gene2 <- k_4_res[rank_num,2]
  gene3 <- k_4_res[rank_num,3]
  gene4 <- k_4_res[rank_num,4]
  # print(paste("current gene1 :", gene1))
  current_score <- matrix(data = 0,ncol = 6,nrow = 2000)
  print(paste("  current time  ",Sys.time(),"current k5 rank  ", rank_num))
  # print(paste("current gene :", gene1,"  -  ",gene2," -  ",gene3))
  for (gene5 in 1:noGene) {
    if(gene5!=gene1&&gene5!=gene2&&gene5!=gene3&&gene5!=gene4){
      score <- sum(pairsZScore[c(gene1,gene2,gene3,gene4,gene5),c(gene1,gene2,gene3,gene4,gene5)])
      if(score > min(current_score[,6])){
        current_score[which.min(current_score[,6]),1:5] = sort(c(gene1,gene2,gene3,gene4,gene5))
        current_score[which.min(current_score[,6]),6] = score
      }
    }
  }
  return(current_score)
}

delOverap <- function(f1,f2){
  item <- c()
  k = dim(f1)[2]
  for (i in 1:dim(f1)[1]) {
    f1Item <- f1[i,]
    for (j in 1:200) {
      a <- intersect(f1Item,f2[j,])
      if(length(a) == k){
        item <- append(item,i)
        break
      }
    }
  }
  if(length(item) != 0){
    f1 <- f1[-item,]
  }
  return(f1)
}

mergeModul <- function(res,cut){
  out_file = sprintf("%s/%s_mergeModul_record.txt", pathOutput,pValue)
  k = dim(res)[2]
  new_modul =  c()
  num_list  = c()
  for(i in 1:(dim(res)[1]-1)){
    if(! i %in% num_list){
      print(paste("current modul: ",i))
      current_modul = res[i,]
      for(j in (i+1):dim(res)[1]){
        next_modul <- res[j,]
        overlap = length(intersect(current_modul,next_modul))
        if(overlap!=0){
          overlap_num = 2*overlap/(length(current_modul)+length(next_modul))
          if(overlap_num>=cut){
            overlap_modul = union(current_modul,next_modul)
            print(as.vector(unlist(current_modul)))
            print(as.vector(unlist(overlap_modul)))
            
            current_score = sum(pairsZScore[as.vector(unlist(current_modul)),as.vector(unlist(current_modul))])/choose(length(current_modul),2)
            overlap_score = sum(pairsZScore[as.vector(unlist(overlap_modul)),as.vector(unlist(overlap_modul))])/choose(length(overlap_modul),2)
            
            x = pairsZScore[as.vector(unlist(current_modul)),as.vector(unlist(current_modul))]
            x_up = x[upper.tri(x)]
            x_sd = sd(x_up)
            # print(paste(" current score: ", current_score))
            # print(paste("overlap score: ",overlap_score))
            current_sd = current_score - x_sd
            cat("current modul: ",as.vector(unlist(current_modul)),"   current score: ", current_score,"   current sd : ",x_sd,"   ",current_sd,"\n",file = out_file,append = T )
            cat("overlap modul: ",as.vector(unlist(overlap_modul)),"   overlap score: ", overlap_score,"\n\n",file = out_file,append = T)
            if(overlap_score>current_sd){
              print(i)
              current_modul = as.vector(unlist(overlap_modul))
              num_list <- append(num_list,j)
            }
          }
        }
      }
      new_modul <- append(new_modul,list(current_modul))
    }
  }
  out = sprintf("%s/%s_mergeModul_%s.txt", pathOutput,pValue,k)
  
  for(i in 1:length(new_modul)){
    print(i)
    cat(as.vector(unlist(new_modul[[i]])),"\n",file = out,sep = '\t',append = T )
  }
}

pathOutput <- "DMMO/output/"
pValue <- "p10"
mutInput <- "DMMO/input/mut.txt"
expInput <- "DMMO/input/exp_pca_10.csv"
oExpInput <- "DMMO/input/RNA_exp.csv"
regInput <- "DMMO/input/node_coReg.txt"
ppiInput <- "DMMO/input/ppi_gae.csv"

#######Calculate mutual exclusion#########
mut <- read.table(mutInput,sep = "\t")
mut <- as.matrix(mut)

noGene <- dim(mut)[2]
mut_colnames <- colnames(mut)
write.table(mut_colnames, file = sprintf("%s/mut_colnames.txt", pathOutput),sep = "\t",quote=FALSE,  row.names = F, col.names = F)


detectCores(logical = F)  
mc <- getOption("mc.cores",80)
res <- mclapply(1:(noGene-1), mut_function, mc.cores = mc)
pairsME <- do.call(rbind, res)
c <- matrix(data = 0,nrow = 1,ncol = noGene)
pairsME <- rbind(pairsME,c)
write.table(pairsME,sprintf("%s/pairsME.txt", pathOutput),sep = "\t")


colnames(pairsME) <- colnames(mut)
rownames(pairsME) <- colnames(mut)

# ME_up = pairsME[upper.tri(pairsME)]
# quantile(ME_up)
# mid = quantile(ME_up)[2]
# pairsME[pairsME < mid] <- 0


zPairsME <- (exp(pairsME^(1/3))-1)/(exp(pairsME^(1/3))+1)
zPairsME <- round(zPairsME,6)
# zPairsME[mutGene,mutGene] <- 0

mut_up <- zPairsME[upper.tri(zPairsME)]
quantile(mut_up)

mid = quantile(mut_up)[3]
zPairsME[zPairsME < mid] <- 0
mut_colnames <-  read.table(sprintf("%s/mut_colnames.txt", pathOutput),sep = "\t",stringsAsFactors=F, check.names=FALSE)


#######Computational co-expression#########

exp_data <- read.table(expInput,sep = ',')
rownames(exp_data) <- exp_data[,1]
exp_data <- exp_data[,-1]
exp_data <- as.matrix(exp_data)
exp <- t(exp_data)

pairsCoExp <- abs(cor(exp,method = "pearson" ))
# pairsCoExp[is.na(pairsCoExp)] <- 0
pairsCoExp[lower.tri(pairsCoExp)] <- 0
diag(pairsCoExp)<-0
pairsCoExp <- round(pairsCoExp,6)

write.table(pairsCoExp,sprintf("%s/pairsCoExp.txt", pathOutput),sep = "\t")


oExp = read.table(oExpInput,sep = ",")

rownames(oExp) = oExp[,1]
oExp = oExp[,-1]
oExp <- t(oExp)

expMean <- colMeans(oExp)
quantile(expMean)

expGene <- colnames(oExp)[which(expMean <  quantile(expMean)[3])]

zPairsCoExp <- pairsCoExp
zPairsCoExp[expGene,expGene] <- 0
# sum(colSums(exp)==0)


######reg#######
REGmatrix<-read.table(regInput,stringsAsFactors = FALSE,sep = "\t")
dim(REGmatrix)

allCoReg <- REGmatrix[,3]
allCoReg <- normalize(allCoReg)
REGmatrix[,3] <- allCoReg
# REGmatrix <- as.matrix(REGmatrix)
write.table(REGmatrix,sprintf("%s/REGMatrix.csv",pathOutput),sep = ',')

pairsCoReg <- matrix(data = 0,ncol = noGene,nrow = noGene)
colnames(pairsCoReg)<-mut_colnames[,1]
rownames(pairsCoReg)<-mut_colnames[,1]

x=1
while (x <= dim(REGmatrix)[1]){
  print(paste("current num:  ",x))
  a<-REGmatrix[x,]
  i<-a[1,1]
  j<-a[1,2]
  k <- a[1,3]
  pairsCoReg[i,j] <- k
  pairsCoReg[j,i] <- k
  x<- x+1
}
write.table(pairsCoReg,sprintf("%s/pairsCoReg.txt", pathOutput),sep = "\t")
pairsCoReg <- round(pairsCoReg,6)


#####ppi#####
ppi <- read.table(ppiInput,sep = ',')
ppi <- t(ppi)
colnames(ppi) <- mut_colnames[,1]
pairsPPI <- abs(cor(ppi,method = "pearson"))
pairsPPI[lower.tri(pairsPPI)] <- 0
pairsPPI <- round(pairsPPI,6)
diag(pairsPPI) <- 0
ppi_up = pairsPPI[upper.tri(pairsPPI)]
quantile(ppi_up)

zPairsPPI <- pairsPPI

pairsZScore<-(pairsCoReg+zPairsME+zPairsPPI+zPairsCoExp)/4
diag(pairsZScore)<-0
pairsZScore = as.matrix(pairsZScore)

ind<-lower.tri(pairsZScore)
pairsZScore[ind] <- 0

write.table(pairsZScore,sprintf("%s/%s_pairsZScore.txt",pathOutput,pValue),sep = "\t")


#k2
sort_matrix <- arrayInd(sort.list(pairsZScore,decreasing=T)[1:10000],dim(pairsZScore))
sort_matrix <- as.matrix(sort_matrix)
write.table(file = sprintf("%s/%s_sort2_index.txt", pathOutput,pValue),x = sort_matrix,sep = '\t',quote = FALSE,row.names = T,col.names = T)


out <- matrix(NA,nrow = 10000,ncol = 2)

out[,1] <- mut_colnames[,1][sort_matrix[,1]]
out[,2] <- mut_colnames[,1][sort_matrix[,2]]

write.table(file = sprintf("%s/%s_sort_2_name.txt",pathOutput,pValue),x = out,sep = '\t',quote = FALSE,row.names = F,col.names = F)

k2Res <- out[1:200,1:2]


print(" k2 finish! ")


k_2_res <- sort_matrix

# mc <- getOption("mc.cores", 50)
res <- mclapply(1:10000, k_3 , mc.cores = mc)
# res <- lapply(1:10000, k_3)
K_3_Score <- do.call(rbind, res)
K_3_Score <- K_3_Score[!duplicated(K_3_Score),]
K_3_Score <- K_3_Score[order(K_3_Score[,4],decreasing = T),][1:10000,]

write.table(file = sprintf("%s/%s_sort_3_index.txt", pathOutput,pValue),x = K_3_Score,sep = '\t',quote = FALSE,row.names = T,col.names = T)


out <- matrix(NA,nrow = 10000,ncol = 4)

out[,1] <- mut_colnames[,1][K_3_Score[,1]]
out[,2] <- mut_colnames[,1][K_3_Score[,2]]
out[,3] <- mut_colnames[,1][K_3_Score[,3]]
out[,4] <- K_3_Score[,4]

k3Res <- out[1:200,1:3]

write.table(file = sprintf("%s/%s_sort_3_name.txt", pathOutput,pValue),x = out,sep = '\t',quote = FALSE,row.names = F,col.names = F)

print("K3 finish")

##k4
k_3_res <- K_3_Score

res <- mclapply(1:2000, k_4 , mc.cores = mc)
# res <- lapply(1:2000, k_4)
K_4_Score <- do.call(rbind, res)
K_4_Score <- K_4_Score[!duplicated(K_4_Score),]
K_4_Score <- K_4_Score[order(K_4_Score[,5],decreasing = T),][1:2000,]


write.table(file = sprintf("%s/%s_sort_4_index.txt", pathOutput,pValue),x = K_4_Score,sep = '\t',quote = FALSE,row.names = T,col.names = T)


out <- matrix(NA,nrow = 2000,ncol = 5)

out[,1] <- mut_colnames[,1][K_4_Score[,1]]
out[,2] <- mut_colnames[,1][K_4_Score[,2]]
out[,3] <- mut_colnames[,1][K_4_Score[,3]]
out[,4] <- mut_colnames[,1][K_4_Score[,4]]
out[,5] <- K_4_Score[,5]

write.table(file = sprintf("%s/%s_sort_4_name.txt", pathOutput,pValue),x = out,sep = '\t',quote = FALSE,row.names = F,col.names = F)

k4Res <- out[1:200,1:4]

print("k4 finish !")



k_4_res <- K_4_Score

# mc <- getOption("mc.cores", 50)
res <- mclapply(1:2000, k_5 , mc.cores = mc)
# res <- lapply(1:2000, k_5)
K_5_Score <- do.call(rbind, res)
K_5_Score <- K_5_Score[!duplicated(K_5_Score),]
K_5_Score <- K_5_Score[order(K_5_Score[,6],decreasing = T),][1:2000,]


write.table(file = sprintf("%s/%s_sort_5_index.txt", pathOutput,pValue),x = K_5_Score,sep = '\t',quote = FALSE,row.names = T,col.names = T)


out <- matrix(NA,nrow = 2000,ncol = 6)

out[,1] <- mut_colnames[,1][K_5_Score[,1]]
out[,2] <- mut_colnames[,1][K_5_Score[,2]]
out[,3] <- mut_colnames[,1][K_5_Score[,3]]
out[,4] <- mut_colnames[,1][K_5_Score[,4]]
out[,5] <- mut_colnames[,1][K_5_Score[,5]]
out[,6] <- K_5_Score[,6]

write.table(file = sprintf("%s/%s_sort_5_name.txt", pathOutput,pValue),x = out,sep = '\t',quote = FALSE,row.names = F,col.names = F)


k5Res <- out[1:200,1:5]

print("k5 finish !")



# del overlap
k2Res <- delOverap(k2Res,k3Res)
k2Res <- delOverap(k2Res,k4Res)
k2Res <- delOverap(k2Res,k5Res)
k3Res <- delOverap(k3Res,k4Res)
k3Res <- delOverap(k3Res,k5Res)
k4Res <- delOverap(k4Res,k5Res)

write.table(k2Res,sprintf("%s/%s_union_name(%s).txt", pathOutput,pValue,2),sep = '\t',quote = FALSE,row.names = F,col.names = F)
write.table(k3Res,sprintf("%s/%s_union_name(%s).txt", pathOutput,pValue,3),sep = '\t',quote = FALSE,row.names = F,col.names = F)
write.table(k4Res,sprintf("%s/%s_union_name(%s).txt", pathOutput,pValue,4),sep = '\t',quote = FALSE,row.names = F,col.names = F)
write.table(k5Res,sprintf("%s/%s_union_name(%s).txt", pathOutput,pValue,5),sep = '\t',quote = FALSE,row.names = F,col.names = F)

print("overlap finish !")


# merge modul
mergeModul(k3Res,0.6)
mergeModul(k4Res,0.7)
mergeModul(k5Res,0.7)

k2Res <- as.data.frame(k2Res)
k3Merge <- read.table(sprintf("%s/%s_mergeModul_3.txt", pathOutput,pValue),fill = T,header = F,sep = "\t")
k4Merge <- read.table(sprintf("%s/%s_mergeModul_4.txt", pathOutput,pValue),fill = T,header = F,sep = "\t")
k5Merge <- read.table(sprintf("%s/%s_mergeModul_5.txt", pathOutput,pValue),fill = T,header = F,sep = "\t")

# com <- rbind.fill(k2Res,k3Merge,k4Merge)
com <- rbind.fill(k2Res,k3Merge,k4Merge,k5Merge)
com[is.na(com)] <- ""
com <- com[!duplicated(com),]

write.table(com,sprintf("%s/com_%s.txt", pathOutput,pValue),sep = '\t',quote = FALSE,row.names = F,col.names = F)

score <- c()

for(i in 1:dim(com)[1]){
  current_case <- c()
  a = com[i,]
  for(j in a){
    if(j!="")
      current_case <- append(current_case,j)
  }
  score[i] <- sum(pairsZScore[current_case,current_case])
  # print(i)
}

write.table(score,sprintf("%s/score_%s.txt", pathOutput,pValue),sep = '\t',quote = FALSE,row.names = F,col.names = F)

print("ALL OVER !")