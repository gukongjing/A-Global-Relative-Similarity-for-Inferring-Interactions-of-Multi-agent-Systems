rm(list = ls())
library(parallel)
#Calculate the number of cores检查电脑当前可用核数
no_cores<-detectCores() - 2		#F-物理CPU核心数/T-逻辑CPU核心数 logical=F

#Initiate cluster发起集群,同时创建数个R进行并行计算
#只是创建待用的核,而不是并行运算环境
cl<-makeCluster(no_cores)

l = c(5,seq(10,310,50))
#现只需要使用并行化版本的lapply,parLapply就可以
final = parLapply(cl, l,function(l){
  rept = 50
  
  traj = function(k,td,Nt,N,uav,A,x0,v0,R){
    library(ggplot2)
    library(network)
    library(GGally)
    library(sna)
    library(RColorBrewer)
    # td =1 
    t <- 0.1
    gnum <- Nt/k
    # R <- 1
    td = td     # time delay
    sdg = 0.01    # 高斯噪声的方差   
    
    # known & unknown parameters
    f = 2
    alpha = 2
    beta = 1
    sigma = 0.5
    
    c <- rep(0,Nt*(N-1)*2)
    w <- rep(0,Nt*(N-1)*2)
    u <- rep(0,Nt*(N-1)*2)
    
    dim(c) <- c(Nt,N-1,2)
    dim(w) <- c(Nt,N-1,2)
    dim(u) <- c(Nt,N-1,2)
    
    
    # initial values
    x <- rep(0,Nt*N*2)
    v <- rep(0,Nt*N*2)
    xg <- rep(0,Nt*N*2)    # 加高斯
    vg <- rep(0,Nt*N*2)    # 加高斯
    
    dim(x) <- c(Nt,N,2)
    dim(v) <- c(Nt,N,2)
    dim(xg) <- c(Nt,N,2)    # 高斯版
    dim(vg) <- c(Nt,N,2)    # 高斯版
    
    x[,(1:(1+td)),] <- x0
    v[,(1:(1+td)),] <- v0
    
    # model
    for(m in (1+td):(N-1)){
      for(j in 1:Nt){
        temp0 <- (x[,m-td,1]-x[j,m,1])^2+(x[,m-td,2]-x[j,m,2])^2
        temp1 <- 1+temp0
        temp2 <- 1/(temp1^f)*A[j,]
        temp31 <- v[,m-td,1]-v[j,m,1]
        temp32 <- v[,m-td,2]-v[j,m,2]
        temp3 <- cbind(temp31,temp32)
        c[j,m,] <- temp2%*%temp3
        
        xm <- colMeans(x[,m-td,],1)
        temp4 <- xm-x[j,m,]
        temp4_ab <- sqrt(sum(temp4^2))
        w[j,m,] <- (1-R/temp4_ab)*temp4
        
        for(i in 1:Nt){
          temp5 <- x[j,m,]-x[i,m-td,]
          if(sqrt(sum(temp5^2)) < 2*sin(pi/Nt)&&sqrt(sum(temp5^2)) != 0){
            temp5_ab <- sqrt(sum(temp5^2))
            temp6 <- temp5/temp5_ab-temp5
            u[j,m,] <- temp6*sin(pi/Nt)/2+u[j,m,]
          }
        }
      }
      
      v[,m+1,] <- (alpha*c[,m,]+beta*u[,m,]+sigma*w[,m,])*t+v[,m,]   #
      x[,m+1,] <- x[,m,]+(v[,m+1,]+v[,m,])*t/2
      
    }
    
    # layout(matrix(c(1,2,3,4),nr=2,byrow = TRUE))
    # layout(1)
    # 
    # plot(x[1,25:N,1],x[1,25:N,2],xlim=c(0,22),ylim=c(0,25),xlab = "x-coordinate",ylab = "y-coordinate",type = 'l',
    #      cex.lab=2,cex.axis=2)   # cex.main=1.5,
    # for (z in 2:Nt) {
    #   lines(x[z,25:N,1],x[z,25:N,2])
    # }
    
    # 添加高斯噪声
    i=0
    for (i in 1:N) {
      temp = sapply(uav, function(uav){
        x[uav,i,] + rnorm(2, mean = 0, sd = sdg)
      })
      xg[,i,] = t(temp)
    }
    
    i=0
    for (i in 1:N) {
      temp = sapply(uav, function(uav){
        v[uav,i,] + rnorm(2, mean = 0, sd = sdg)
      })
      vg[,i,] = t(temp)
    }
    
    resu= list(xg,vg)
    return(resu)
  }
  
  
  
  ####similarity measurement
  library(dtwclust)
  library(dtw)
  library(forcats)
  library(igraph)
  library(SimilarityMeasures)
  library(infotheo)
  library(ggm)
  library(energy)
  
  #### 2 固定k=5
  Kfixed5 = function(w,series){    #k=clustering number, w=historical impact
    
    kk = 5      #k固定为5
    # store the clustering result
    division = matrix(0, nrow = length(sequence), ncol = Nt)                     # result of each piece
    C = matrix(0,nrow = Nt, ncol = Nt) 
    i=0
    for(i in sequence){
      piece = t(sapply(uav,function(uav){as.matrix(series[[uav]][i,])}))      # trajectory
      clus = as.numeric(fct_inorder(as.factor(kmeans(piece,centers = kk)$cluster)))
      
      D = matrix(0,nrow = Nt, ncol = Nt)             # initialize before each for
      for (ki in 1:kk) {
        index = which(clus == ki)
        D[index,index] = 1           # count for each piece of network
      }
      C = D + w*C         # accumulated network
      
      division[i,] = clus
    }
    
    simimatrix = C
    diag(simimatrix) = 0
    simimatrix = simimatrix/(max(simimatrix))
    return(simimatrix)
  }
  
  #### 3 固定k=8
  Kfixed8 = function(w,series){    #k=clustering number, w=historical impact
    
    kk = 8      #k固定为8
    # store the clustering result
    division = matrix(0, nrow = length(sequence), ncol = Nt)                     # result of each piece
    C = matrix(0,nrow = Nt, ncol = Nt) 
    i=0
    for(i in sequence){
      piece = t(sapply(uav,function(uav){as.matrix(series[[uav]][i,])}))      # trajectory
      clus = as.numeric(fct_inorder(as.factor(kmeans(piece,centers = kk)$cluster)))
      
      D = matrix(0,nrow = Nt, ncol = Nt)             # initialize before each for
      for (ki in 1:kk) {
        index = which(clus == ki)
        D[index,index] = 1           # count for each piece of network
      }
      C = D + w*C         # accumulated network
      
      division[i,] = clus
    }
    
    simimatrix = C
    diag(simimatrix) = 0
    simimatrix = simimatrix/(max(simimatrix))
    return(simimatrix)
  }
  
  ### 4 随机K   多变量
  Krandom = function(w,series){    #k=clustering number, w=historical impact
    
    # store the clustering result
    division = matrix(0, nrow = length(sequence), ncol = Nt)                     # result of each piece
    C = matrix(0,nrow = Nt, ncol = Nt) 
    i=1
    for(i in sequence){
      kk = sample(2:(Nt-1),1)
      piece = t(sapply(uav,function(uav){as.matrix(series[[uav]][i,])}))      # trajectory
      clus = as.numeric(fct_inorder(as.factor(kmeans(piece,centers = kk)$cluster)))
      
      D = matrix(0,nrow = Nt, ncol = Nt)             # initialize before each for
      for (ki in 1:kk) {
        index = which(clus == ki)
        D[index,index] = 1           # count for each piece of network
      }
      C = D + w*C         # accumulated network
      
      division[i,] = clus
    }
    
    simimatrix = C
    diag(simimatrix) = 0
    simimatrix = simimatrix/max(simimatrix)
    return(simimatrix)
  }
  
  ### 5 Change detection
  krandomNN = function(w,series){    #k=clustering number, w=historical impact,nd=net distance
    # store the clustering result
    C = matrix(0,Nt,Nt)
    i=0
    for(i in sequence){
      piece = t(sapply(uav,function(uav){as.matrix(series[[uav]][i,])}))      # trajectory
      D = matrix(0,Nt,Nt)
      for (j in 1:floor(Nt/2)) {
        kk = sample(2:(Nt-1),1)
        E = matrix(0, Nt, Nt)
        clus = as.numeric(fct_inorder(as.factor(kmeans(piece,centers = kk)$cluster)))
        for (ki in 1:kk) {
          index = which(clus == ki)
          E[index,index] = 1           # count for each piece of network
        }
        diag(E) = 0
        D = D + E
      }
      C = D + w*C
    }
    simimatrix = C
    diag(simimatrix) = 0
    simimatrix = simimatrix/(max(simimatrix))
    return(simimatrix)
  }
  
  
  ### 7 Change detection
  KfixedN = function(w,series){    #k=clustering number, w=historical impact,nd=net distance
    # store the clustering result
    C = matrix(0,Nt,Nt)
    kk = floor(Nt/2)
    i=0
    for(i in sequence){
      piece = t(sapply(uav,function(uav){as.matrix(series[[uav]][i,])}))      # trajectory
      D = matrix(0,Nt,Nt)
      clus = as.numeric(fct_inorder(as.factor(kmeans(piece,centers = kk)$cluster)))
      for (ki in 1:kk) {
        index = which(clus == ki)
        D[index,index] = 1           # count for each piece of network
      }
      diag(D) = 0
      C = C+w*D
    }
    
    simimatrix = C
    diag(simimatrix) = 0
    simimatrix = simimatrix/(max(simimatrix))
    return(simimatrix)
  }
  
  #### 8
  Kadaptive = function(w,series){    #k=clustering number, w=historical impact
    
    # store the clustering result
    division = matrix(0, nrow = length(sequence), ncol = Nt)                     # result of each piece
    C = matrix(0,nrow = Nt, ncol = Nt) 
    i=0
    for(i in sequence){
      piece = t(sapply(uav,function(uav){as.matrix(series[[uav]][i,])}))      # trajectory
      nk = c(2:10)
      SW =sapply(nk, function(nk){                                           #计算使用kmeans方法时每个待测试的nk的轮廓系数
        cluster.stats(dist(piece),kmeans(piece,centers = nk, nstart = 3)$cluster)$avg.silwidth
      })
      kk = nk[which.max(SW)]   
      
      clus = as.numeric(fct_inorder(as.factor(kmeans(piece,centers = kk)$cluster)))
      
      D = matrix(0,nrow = Nt, ncol = Nt)             # initialize before each for
      for (ki in 1:kk) {
        index = which(clus == ki)
        D[index,index] = 1           # count for each piece of network
      }
      diag(D) = 0
      C = D + w*C         # accumulated network
    }
    
    simimatrix = C
    diag(simimatrix) = 0
    simimatrix = simimatrix/(max(simimatrix))
    return(simimatrix)
  }
  
  
  library(igraph)
  library(network)
  library(fpc)
  # library(minet)
  
  
  prauc = matrix(0, nrow = rept, ncol = 6)
  rocauc = matrix(0, nrow = rept, ncol = 6)
  
  for (rr in 1:rept) {
    
    Nt <- 40   #数量
    uav = c(1:Nt)   #无人机索引
    kind = 3        # 网络类型
    k1 = 5           #子群数目
    k2 = 8
    m = 3
    disturbnum = 400
    td = 1
    N = 600
    
    # 生成参考网络
    ## 基础矩阵:类对角矩阵
    temp = Nt-1
    Atemp = diag(temp)
    A = matrix(0,nrow = Nt, ncol = Nt)     #类对角
    A[2:Nt,1:(Nt-1)] = Atemp
    A[upper.tri(A)] = t(A)[upper.tri(t(A))]
    LDJ = A     #类对角
    ## ER network
    p = 0.01
    ernet = erdos.renyi.game(Nt, p, type=c("gnp"), directed = FALSE, loops = FALSE)
    # plot(ernet)
    ernet = as.matrix(get.adjacency(ernet))
    gnet = ernet + LDJ
    gnet[which(gnet>1)] = 1
    diag(gnet)=0
    
    # G = graph.adjacency(gnet,mode="undirected")
    # plot(G)
    
    # 生成轨迹信息
    x0 <- runif(Nt*2*(1+td),min = 0,max = 10)
    v0 <- runif(Nt*2*(1+td),min = 0,max = 1)
    dim(x0) = c(Nt,(td+1),2)
    dim(v0) = c(Nt,(1+td),2)
    xvg = traj(k=2,td,Nt,N,uav,gnet,x0,v0,R=1)
    
    xg = xvg[[1]]
    vg = xvg[[2]]
    
    plot(xg[Nt,,1],xg[Nt,,2],type = "l") #,xlim=c(-20,20),ylim=c(-20,20)
    for (i0 in 1:(Nt-1)) {
      lines(xg[i0,,1],xg[i0,,2])
    }
    
    sequence = c(1:l)
    
    start = l+1
    sam = seq(start,N-1,floor(l/3*2))
    auc = rep(0,length(sam)*6*2)
    dim(auc) = c(length(sam),6,2)
    i = 0
    
    for (time in sam) {
      i = i+1
      temp = c((time-l):(time-1))
      
      xstand1 = (xg[,temp,1]-mean(xg[,temp,1]))/sqrt(sum((xg[,temp,1]-mean(xg[,temp,1]))^2)/(length(xg[,temp,1])-1))
      vstand1 = (vg[,temp,1]-mean(vg[,temp,1]))/sqrt(sum((vg[,temp,1]-mean(vg[,temp,1]))^2)/(length(vg[,temp,1])-1))
      xstand2 = (xg[,temp,2]-mean(xg[,temp,2]))/sqrt(sum((xg[,temp,2]-mean(xg[,temp,2]))^2)/(length(xg[,temp,2])-1))
      vstand2 = (vg[,temp,2]-mean(vg[,temp,2]))/sqrt(sum((vg[,temp,2]-mean(vg[,temp,2]))^2)/(length(vg[,temp,2])-1))
      
      
      series1xv = lapply(uav, function(uav){
        as.matrix(data.frame(coordsx1 = xstand1[uav,], coordsx2 = xstand2[uav,],
                             coordsv1 = vstand1[uav,], coordsv2 = vstand2[uav,]))
      })
      series = series1xv
      
      
      library(ROCR)
      library(ROCR)
      
      w=1
      library(corrplot)
      ### similarity
      
      #2 
      Kfixed5simi = Kfixed5(w=1,series)

      pred = prediction(as.vector(Kfixed5simi[upper.tri(Kfixed5simi)]),as.vector(gnet[upper.tri(gnet)]))
      pr.perf = performance(pred, measure = "aucpr")
      auc[i,1,1] = unlist(pr.perf@y.values)
      auc.perf = performance(pred, measure = "auc")
      auc[i,1,2] = unlist(auc.perf@y.values)
      
      #3 
      Kfixed8simi = Kfixed8(w,series)

      pred = prediction(as.vector(Kfixed8simi[upper.tri(Kfixed8simi)]),as.vector(gnet[upper.tri(gnet)]))
      pr.perf = performance(pred, measure = "aucpr")
      auc[i,2,1] = unlist(pr.perf@y.values)
      auc.perf = performance(pred, measure = "auc")
      auc[i,2,2] = unlist(auc.perf@y.values)
      
      #4 
      Krandomsimi = Krandom(w,series)

      pred = prediction(as.vector(Krandomsimi[upper.tri(Krandomsimi)]),as.vector(gnet[upper.tri(gnet)]))
      pr.perf = performance(pred, measure = "aucpr")
      auc[i,3,1] = unlist(pr.perf@y.values)
      auc.perf = performance(pred, measure = "auc")
      auc[i,3,2] = unlist(auc.perf@y.values)
      
      #5
      kmultisimi = krandomNN(w,series)
 
      pred = prediction(as.vector(kmultisimi[upper.tri(kmultisimi)]),as.vector(gnet[upper.tri(gnet)]))
      pr.perf = performance(pred, measure = "aucpr")
      auc[i,4,1] = unlist(pr.perf@y.values)
      auc.perf = performance(pred, measure = "auc")
      auc[i,4,2] = unlist(auc.perf@y.values)
      
      #7
      KmultiNsimi = KfixedN(w,series)

      pred = prediction(as.vector(KmultiNsimi[upper.tri(KmultiNsimi)]),as.vector(gnet[upper.tri(gnet)]))
      pr.perf = performance(pred, measure = "aucpr")
      auc[i,5,1] = unlist(pr.perf@y.values)
      auc.perf = performance(pred, measure = "auc")
      auc[i,5,2] = unlist(auc.perf@y.values)
      
      #8
      Kadaptivesimi = Kadaptive(w,series)

      pred = prediction(as.vector(Kadaptivesimi[upper.tri(Kadaptivesimi)]),as.vector(gnet[upper.tri(gnet)]))
      pr.perf = performance(pred, measure = "aucpr")
      auc[i,6,1] = unlist(pr.perf@y.values)
      auc.perf = performance(pred, measure = "auc")
      auc[i,6,2] = unlist(auc.perf@y.values)
      
    }
    
    prauc[rr,] = apply(auc[,,1], 2, mean)
    rocauc[rr,] = apply(auc[,,2], 2, mean)
    
  }
  
  D = c(apply(prauc, 2, mean),apply(rocauc, 2, mean))
  return(D)
})

stopCluster(cl)

D = do.call(rbind, final)

# write.table(D,file = 'different_k_windowsize.txt', sep = ' ', row.names = FALSE, col.names = FALSE)
D = as.matrix(read.table('different_k_windowsize.txt',sep = ' ', header = FALSE))

E = as.vector(D)
b = length(E)
prauc = E[1:(b/2)]
rocauc = E[(b/2+1):b]

prmean = matrix(prauc,nrow = length(l))
prmean = round(apply(prmean, 2, mean),2)
rocmean = matrix(rocauc,nrow = length(l))
rocmean = round(apply(rocmean, 2, mean),2)

library(ggplot2)
library(RColorBrewer)
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
color1 <- colorRampPalette(rev(brewer.pal(11,'Spectral')))(28)
color2 = c("#dc745c","#f2c368",'#bdd7ba',"#b0bcd6","#6a84b5","#3b5899",'#041656')  #"#efaa6e",,

###################################四图并列###############################################################

ggdata11 = data.frame(PRAUC = prmean, Method = c('K-fixed=5','K-fixed=8','K-random','K-random*N/2',
                                               'K-fixed=N/2','K-adaptive'))
ggdata11$Method = factor(ggdata11$Method,levels = c('K-random*N/2','K-random',
                                                    'K-adaptive','K-fixed=N/2','K-fixed=8','K-fixed=5'))
p11 = ggplot(data = ggdata11, mapping = aes(x = Method, y =PRAUC))+
  theme(axis.title.x = element_text(size = 12,family='sans'), axis.text.x = element_text(angle = 35, hjust=1, size = 12,family='sans'),
        axis.title.y = element_text(size = 12,family='sans'), axis.text.y = element_text(size = 12,family='sans'),
        title=element_text(size=15,family='sans'),legend.text = element_text(size = 12,family='sans'))+ 
  geom_bar(stat="identity",aes(fill = Method))+guides(fill='none')+ coord_cartesian(ylim = c(0.25, 0.6))+
  scale_fill_manual(values=color2)+ylab("PR AUC")
p11
ggdata12 = data.frame(PRAUC = prauc, Method = rep(c('K-fixed=5','K-fixed=8','K-random','K-random*N/2',
                                                  'K-fixed=N/2','K-adaptive'),each = length(l)),
                      Windowsize = rep(l,6))
ggdata12$Method = factor(ggdata12$Method,levels = c('K-random*N/2','K-random',
                                                    'K-adaptive','K-fixed=N/2','K-fixed=8','K-fixed=5'))

p12 = ggplot(data = ggdata12, mapping = aes(x = Windowsize, y = PRAUC,color = Method))+
  geom_point() + geom_line()+
  theme(axis.title.x = element_text(size = 12,family = 'sans'), axis.text.x = element_text(angle = 45, hjust=1, 
                                                                                           size = 12,family = 'sans'),
        axis.title.y = element_text(size = 12,family = 'sans'), axis.text.y = element_text(size = 12,family = 'sans'),
        title=element_text(size=12,family = 'sans'),legend.position = 'none')+ 
  scale_color_manual(values=color2)+ylab("PR AUC")
p12
library(cowplot)
gg1 <- ggdraw() + 
  draw_plot(p12, 0,0,0.5, 1) + draw_plot(p11, 0.5,0,0.5,1)
print(gg1)



ggdata21 = data.frame(ROCAUC = rocmean, Method = c('K-fixed=5','K-fixed=8','K-random','K-random*N/2',
                                                'K-fixed=N/2','K-adaptive'))
ggdata21$Method = factor(ggdata21$Method,levels = c('K-random*N/2','K-random',
                                                    'K-adaptive','K-fixed=N/2','K-fixed=8','K-fixed=5'))
p21 = ggplot(data = ggdata21, mapping = aes(x = Method, y = ROCAUC))+
  theme(axis.title.x = element_text(size = 12,family='sans'), axis.text.x = element_text(angle = 35, hjust=1, size = 12,family='sans'),
        axis.title.y = element_text(size = 12,family='sans'), axis.text.y = element_text(size = 12,family='sans'),
        title=element_text(size=15,family='sans'),legend.text = element_text(size = 12,family='sans'))+ 
  geom_bar(stat="identity",aes(fill = Method))+guides(fill='none')+ coord_cartesian(ylim = c(0.8, 0.9))+
  scale_y_continuous(breaks=c(0.8,0.85,0.9))+
  scale_fill_manual(values=color2)+ylab("ROC AUC")
p21
ggdata22 = data.frame(ROCAUC = rocauc, Method = rep(c('K-fixed=5','K-fixed=8','K-random','K-random*N/2',
                                                   'K-fixed=N/2','K-adaptive'),each = length(l)),
                      Windowsize = rep(l,6))
ggdata22$Method = factor(ggdata22$Method,levels = c('K-random*N/2','K-random',
                                                    'K-adaptive','K-fixed=N/2','K-fixed=8','K-fixed=5'))

p22 = ggplot(data = ggdata22, mapping = aes(x = Windowsize, y = ROCAUC,color = Method))+
  geom_point() + geom_line()+
  theme(axis.title.x = element_text(size = 12,family = 'sans'), axis.text.x = element_text(angle = 45, hjust=1, 
                                                                                           size = 12,family = 'sans'),
        axis.title.y = element_text(size = 12,family = 'sans'), axis.text.y = element_text(size = 12,family = 'sans'),
        title=element_text(size=12,family = 'sans'),legend.position = 'none')+ 
  scale_color_manual(values=color2)+ylab("ROC AUC")
p22
library(cowplot)
gg2 <- ggdraw() + 
  draw_plot(p22, 0,0,0.5, 1) + draw_plot(p21, 0.5,0,0.5,1)
print(gg2)

library(cowplot)
gg <- ggdraw() + 
  draw_plot(gg1, 0,0,0.5, 1) + draw_plot(gg2, 0.5,0,0.5,1)
print(gg)

######################################两图并列#####################################################
# ggdata11 = data.frame(AUC = prmean, Method = c('Kfixed=2','Kfixed=5','Kfixed=8','Krandom','Krandom*N/2',
#                                                'Krandom*N','Kfixed=N/2','Kadaptive'))
# ggdata11$Method = factor(ggdata11$Method,levels = c('Krandom*N','Krandom*N/2','Krandom',
#                                                     'Kadaptive','Kfixed=N/2','Kfixed=8','Kfixed=5','Kfixed=2'))
# p11 = ggplot(data = ggdata11, mapping = aes(x = Method, y = AUC,label = round(AUC,3)))+
#   theme(axis.title.x = element_text(size = 12,family='sans'), axis.text.x = element_text(angle = 35, hjust=1, size = 12,family='sans'),
#         axis.title.y = element_text(size = 12,family='sans'), axis.text.y = element_text(size = 12,family='sans'),
#         title=element_text(size=15,family='sans'),legend.text = element_text(size = 12,family='sans'))+ 
#   geom_text(aes(y = AUC), position = position_dodge(0.9), vjust = -0.5)+
#   geom_bar(stat="identity",aes(fill = Method))+guides(fill='none')+ coord_cartesian(ylim = c(0.15, 0.5))+
#   scale_fill_manual(values=color2)
# p11
# ggdata12 = data.frame(AUC = prauc, Method = rep(c('Kfixed=2','Kfixed=5','Kfixed=8','Krandom','Krandom*N/2',
#                                                   'Krandom*N','Kfixed=N/2','Kadaptive'),each = length(l)),
#                       Windowsize = rep(l,8))
# ggdata12$Method = factor(ggdata12$Method,levels = c('Krandom*N','Krandom*N/2','Krandom',
#                                                     'Kadaptive','Kfixed=N/2','Kfixed=8','Kfixed=5','Kfixed=2'))
# 
# p12 = ggplot(data = ggdata12, mapping = aes(x = Windowsize, y = AUC,color = Method))+
#   geom_point() + geom_line()+ggtitle("PR AUC") +
#   theme(axis.title.x = element_text(size = 12,family = 'sans'), axis.text.x = element_text(angle = 45, hjust=1, 
#                                                                                            size = 12,family = 'sans'),
#         axis.title.y = element_text(size = 12,family = 'sans'), axis.text.y = element_text(size = 12,family = 'sans'),
#         title=element_text(size=12,family = 'sans'),legend.position = 'none')+ 
#   coord_cartesian(ylim = c(-0.1, 0.5))+
#   scale_color_manual(values=color2)
# p12
# library(cowplot)
# gg1 <- ggdraw() + 
#   draw_plot(p12, 0,0,1, 1) + draw_plot(p11, 0.3,0.12,0.68,0.4)
# print(gg1)
# 
# 
# 
# ggdata21 = data.frame(AUC = rocmean, Method = c('Kfixed=2','Kfixed=5','Kfixed=8','Krandom','Krandom*N/2',
#                                                 'Krandom*N','Kfixed=N/2','Kadaptive'))
# ggdata21$Method = factor(ggdata21$Method,levels = c('Krandom*N','Krandom*N/2','Krandom',
#                                                     'Kadaptive','Kfixed=N/2','Kfixed=8','Kfixed=5','Kfixed=2'))
# p21 = ggplot(data = ggdata21, mapping = aes(x = Method, y = AUC,label = round(AUC,3)))+
#   theme(axis.title.x = element_text(size = 12,family='sans'), axis.text.x = element_text(angle = 35, hjust=1, size = 12,family='sans'),
#         axis.title.y = element_text(size = 12,family='sans'), axis.text.y = element_text(size = 12,family='sans'),
#         title=element_text(size=15,family='sans'),legend.text = element_text(size = 12,family='sans'))+ 
#   geom_text(aes(y = AUC), position = position_dodge(0.9), vjust = -0.5)+
#   geom_bar(stat="identity",aes(fill = Method))+guides(fill='none')+ coord_cartesian(ylim = c(0.7, 0.82))+
#   scale_y_continuous(breaks=c(0.7,0.75,0.8))+
#   scale_fill_manual(values=color2)
# p21
# ggdata22 = data.frame(AUC = rocauc, Method = rep(c('Kfixed=2','Kfixed=5','Kfixed=8','Krandom','Krandom*N/2',
#                                                    'Krandom*N','Kfixed=N/2','Kadaptive'),each = length(l)),
#                       Windowsize = rep(l,8))
# ggdata22$Method = factor(ggdata22$Method,levels = c('Krandom*N','Krandom*N/2','Krandom',
#                                                     'Kadaptive','Kfixed=N/2','Kfixed=8','Kfixed=5','Kfixed=2'))
# 
# p22 = ggplot(data = ggdata22, mapping = aes(x = Windowsize, y = AUC,color = Method))+
#   geom_point() + geom_line(size= 1,alpha = 0.5)+ggtitle("ROC AUC") +
#   theme(axis.title.x = element_text(size = 12,family = 'sans'), axis.text.x = element_text(angle = 45, hjust=1, 
#                                                                                            size = 12,family = 'sans'),
#         axis.title.y = element_text(size = 12,family = 'sans'), axis.text.y = element_text(size = 12,family = 'sans'),
#         title=element_text(size=12,family = 'sans'),legend.position = 'none')+ 
#   coord_cartesian(ylim = c(0.6, 0.82))+
#   scale_color_manual(values=color2)
# p22
# library(cowplot)
# gg2 <- ggdraw() + 
#   draw_plot(p22, 0,0,1, 1) + draw_plot(p21, 0.3,0.12,0.68,0.4)
# print(gg2)
# 
# library(cowplot)
# gg <- ggdraw() + 
#   draw_plot(gg1, 0,0,0.48, 1) + draw_plot(gg2, 0.52,0,0.48,1)+
#   draw_plot_label(c("A", "B"), c(0, 0.52), c(1,1), size = 15,colour = "black")
# print(gg)