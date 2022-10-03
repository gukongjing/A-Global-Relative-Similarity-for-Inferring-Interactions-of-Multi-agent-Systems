rm(list = ls())
library(parallel)
#Calculate the number of cores检查电脑当前可用核数
no_cores<-detectCores() - 2		#F-物理CPU核心数/T-逻辑CPU核心数 logical=F

#Initiate cluster发起集群,同时创建数个R进行并行计算
#只是创建待用的核,而不是并行运算环境
cl<-makeCluster(no_cores)

rept = c(1:50)
#现只需要使用并行化版本的lapply,parLapply就可以
final = parLapply(cl, rept,function(rept){
  
  communication = 2.5
  
  Nt = 40   #数量
  uav = c(1:Nt)   #无人机索引
  td = sample(1:3, 1, replace = FALSE)
  N = 301
  # set.seed(5*rept+72)
  
  traj = function(k,td,Nt,N,uav,communication,x0,v0){
    library(ggplot2)
    library(network)
    library(GGally)
    library(sna)
    library(RColorBrewer)
    # td =1 
    t <- 0.1
    gnum <- Nt/k
    R <- 3
    td = td     # time delay
    sdg = 0.01    # 高斯噪声的方差
    out = 0.1
    
    
    # known & unknown parameters
    f = 2
    alpha = 1
    beta = 0.5
    sigma = 0.2
    
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
    
    # define conection matrix
    E = array(0,dim = c(Nt,Nt,N-1))
    A = matrix(0,nrow=Nt,ncol=Nt)
    
    # model
    for(m in (1+td):(N-1)){
      
      # setting network connection
      for (no in 1:Nt) {   
        TR = (apply((x[no,m,]-t(x[,m,]))^2,2,sum) < communication^2)     #  quadratic sum,less than 10
        A[no,] = as.numeric(TR)
      } 
      
      diag(A)=0
      E[,,m] = A
      
      
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
            u[j,m,] <- temp6*R*sin(pi/Nt)/2+u[j,m,]
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
    
    
    resu= list(xg,vg,E)
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
  
  ### 1 ED
  EDsimi = function(series){
    #计算距离矩阵
    i=0
    simimatrix = matrix(0,nrow = Nt, ncol = Nt)
    for (i in 1:Nt) {
      j=0
      for (j in 1:Nt) {
        simimatrix[i,j] =
          mean(sapply(sequence, function(sequence){
            sqrt(sum((series[[i]][sequence,] - series[[j]][sequence,])^2))
          }))
      }
    }
    
    simimatrix = max(simimatrix)-simimatrix
    diag(simimatrix) = 0
    simimatrix = simimatrix/max(simimatrix)
    return(simimatrix)
  }
  
  
  ### 2 DTW
  DTWsimi = function(series){
    
    #计算距离矩阵
    i=0
    simimatrix = matrix(0,nrow = Nt, ncol = Nt)
    for (i in 1:Nt) {
      simimatrix[i,] = sapply(uav, function(uav){
        # Multivariate series must have time spanning the rows and variables spanning the columns.
        # 行为时间，列为变量
        dtw_basic(series[[i]], y = series[[uav]], window.size = NULL, norm = "L2")
      })
    }
    
    simimatrix = max(simimatrix)-simimatrix
    diag(simimatrix) = 0
    simimatrix = simimatrix/max(simimatrix)
    return(simimatrix)
  }
  
  
  
  ### 3 随机K   多变量
  Krandom = function(w,series){    #k=clustering number, w=historical impact
    
    # store the clustering result
    division = matrix(0, nrow = length(sequence), ncol = Nt)                     # result of each piece
    C = matrix(0,nrow = Nt, ncol = Nt) 
    i=1
    for(i in sequence){
      piece = t(sapply(uav,function(uav){as.matrix(series[[uav]][i,])}))      # trajectory
      # va = mean(sapply(uav, function(i){
      #   mean(sapply(uav,function(j){
      #     sum(sqrt((piece[i,]-piece[j,])^2))
      #   }))
      # }))
      for (j in 1:(Nt)) {
        kk = sample(2:(Nt-1),1)
        clus = as.numeric(fct_inorder(as.factor(kmeans(piece,centers = kk)$cluster)))
        D = matrix(0,nrow = Nt, ncol = Nt)             # initialize before each for
        for (ki in 1:kk) {
          index = which(clus == ki)
          D[index,index] = 1           # count for each piece of network
        }
        C = D + C         # accumulated network
      }
      
    }
    
    simimatrix = C
    diag(simimatrix) = 0
    simimatrix = simimatrix/max(simimatrix)
    return(simimatrix)
  }
  
  
  ### 4 mutual information
  MIsimi = function(series){
    i = 0 
    simimatrix = matrix(0, nrow = Nt, ncol = Nt)
    for (i in uav) {
      simimatrix[i,] = sapply(uav, function(uav){
        ts = discretize(rbind(series[[i]],series[[uav]]))
        mutinformation(ts[1:(dim(ts)[1]/2),], ts[(dim(ts)[1]/2+1):(dim(ts)[1]),], method="emp")
      })
    }
    diag(simimatrix) = 0
    simimatrix = simimatrix/max(simimatrix)
    return(simimatrix)
  }
  
  ### 5 Distance correlation
  DCsimi = function(series){
    i = 0 
    simimatrix = matrix(0, nrow = Nt, ncol = Nt)
    for (i in uav) {
      simimatrix[i,] = sapply(uav,function(uav){dcor(series[[i]],series[[uav]])})
    }
    diag(simimatrix) = 0
    simimatrix = simimatrix/max(simimatrix)
    return(simimatrix)
  }
  
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
  
  
  library(igraph)
  library(network)
  # library(minet)
  
  
  # 生成轨迹信息
  x0 <- runif(Nt*2*(1+td),min = 0,max = 10)
  v0 <- runif(Nt*2*(1+td),min = 0,max = 1)
  dim(x0) = c(Nt,(td+1),2)
  dim(v0) = c(Nt,(1+td),2)
  xvE = traj(k=2,td,Nt,N,uav,communication,x0,v0)
  xg = xvE[[1]]
  vg = xvE[[2]]
  
  plot(xvE[[1]][Nt,,1],xvE[[1]][Nt,,2],type = "l") #,xlim=c(-20,20),ylim=c(-20,20)
  for (i0 in 1:(Nt-1)) {
    lines(xvE[[1]][i0,,1],xvE[[1]][i0,,2])
  }
  
  
  l = 5
  sequence = c(1:(l+1))
  
  start = 6
  identification = seq(start,(N-1),20)
  auc = rep(0,length(identification)*6*2)
  dim(auc) = c(length(identification),6,2)
  i = 0
  
  for (time in identification) {
    i = i+1
    temp = c((time-l):time)
    
    xstand1 = (xg[,temp,1]-mean(xg[,temp,1]))/sqrt(sum((xg[,temp,1]-mean(xg[,temp,1]))^2)/(length(xg[,temp,1])-1))
    vstand1 = (vg[,temp,1]-mean(vg[,temp,1]))/sqrt(sum((vg[,temp,1]-mean(vg[,temp,1]))^2)/(length(vg[,temp,1])-1))
    xstand2 = (xg[,temp,2]-mean(xg[,temp,2]))/sqrt(sum((xg[,temp,2]-mean(xg[,temp,2]))^2)/(length(xg[,temp,2])-1))
    vstand2 = (vg[,temp,2]-mean(vg[,temp,2]))/sqrt(sum((vg[,temp,2]-mean(vg[,temp,2]))^2)/(length(vg[,temp,2])-1))
    
    
    series1xv = lapply(uav, function(uav){
      as.matrix(data.frame(coordsx1 = xstand1[uav,], coordsx2 = xstand2[uav,],
                           coordsv1 = vstand1[uav,], coordsv2 = vstand2[uav,]))
    })
    series = series1xv
    
    gnet = xvE[[3]][,,time]
    
    library(ROCR)
    
    ## 1 ED
    edsimilarity = EDsimi(series)
    # corrplot(edsimilarity)
    
    pred = prediction(c(as.vector(edsimilarity[upper.tri(edsimilarity)]),0),
                      c(as.vector(gnet[upper.tri(gnet)]),0))
    pr.perf = performance(pred, measure = "aucpr")
    auc[i,1,1] = unlist(pr.perf@y.values)
    auc.perf = performance(pred, measure = "auc")
    auc[i,1,2] = unlist(auc.perf@y.values)
    
    #2 dtw
    dtwsimilarity = DTWsimi(series)
    
    pred = prediction(c(as.vector(dtwsimilarity[upper.tri(dtwsimilarity)]),0),
                      c(as.vector(gnet[upper.tri(gnet)]),0))
    pr.perf = performance(pred, measure = "aucpr")
    auc[i,2,1] = unlist(pr.perf@y.values)
    auc.perf = performance(pred, measure = "auc")
    auc[i,2,2] = unlist(auc.perf@y.values)
    
    #3 Ksimi
    w=1  #historical impact
    Ksimilarity = Krandom(w,series)
    
    pred = prediction(c(as.vector(Ksimilarity[upper.tri(Ksimilarity)]),0),
                      c(as.vector(gnet[upper.tri(gnet)]),0))
    pr.perf = performance(pred, measure = "aucpr")
    auc[i,3,1] = unlist(pr.perf@y.values)
    auc.perf = performance(pred, measure = "auc")
    auc[i,3,2] = unlist(auc.perf@y.values)
    
    #4 MI
    MIsimilarity = MIsimi(series)
    
    pred = prediction(c(as.vector(MIsimilarity[upper.tri(MIsimilarity)]),0),
                      c(as.vector(gnet[upper.tri(gnet)]),0))
    pr.perf = performance(pred, measure = "aucpr")
    auc[i,4,1] = unlist(pr.perf@y.values)
    auc.perf = performance(pred, measure = "auc")
    auc[i,4,2] = unlist(auc.perf@y.values)
    
    #5 pcov
    pcorsimilarity = DCsimi(series)
    
    pred = prediction(c(as.vector(pcorsimilarity[upper.tri(pcorsimilarity)]),0),
                      c(as.vector(gnet[upper.tri(gnet)]),0))
    pr.perf = performance(pred, measure = "aucpr")
    auc[i,5,1] = unlist(pr.perf@y.values)
    auc.perf = performance(pred, measure = "auc")
    auc[i,5,2] = unlist(auc.perf@y.values)
    
    #6
    kmultisimi = krandomNN(w,series)
    
    pred = prediction(as.vector(kmultisimi[upper.tri(kmultisimi)]),as.vector(gnet[upper.tri(gnet)]))
    pr.perf = performance(pred, measure = "aucpr")
    auc[i,6,1] = unlist(pr.perf@y.values)
    auc.perf = performance(pred, measure = "auc")
    auc[i,6,2] = unlist(auc.perf@y.values)
  }
  
  prauc = auc[,,1]
  rocauc = auc[,,2]
  
  D = c(as.vector(prauc),as.vector(rocauc))
  
  return(D)
})

stopCluster(cl)

D = do.call(rbind, final)

# write.table(D,file = 'time_varying_comparison.txt', sep = ' ', row.names = FALSE, col.names = FALSE)
D = as.matrix(read.table('time_varying_comparison.txt',sep = ' ', header = FALSE))

E = apply(D, 2, mean)
b = length(E)

prauc = E[1:(b/2)]
rocauc = E[(b/2+1):b]
time = seq(6,300,20)

library(RColorBrewer)
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
color1 <- colorRampPalette(rev(brewer.pal(11,'Spectral')))(28)
# color2 = c("#dc745c","#efaa6e","#f2c368","#6a84b5","#3b5899")   #"#b0bcd6",
color2 = c("#dc745c","#f2c368","#b0bcd6","#6a84b5","#3b5899",'#041656')  #"#efaa6e",

# 画图AUC-Netdens
ggdatadens = data.frame(AUC = c(as.vector(prauc),as.vector(rocauc)), 
                        Method = rep(rep(c('ED','DTW','Krandom','MI','DC','Krandom*N/2'),
                                         each = length(time)),2),
                        Time = rep(rep(c(time), 6),2),
                        Type = rep(c("PR",'ROC'),each = length(time)*6))
ggdatadens$Method <- factor(ggdatadens$Method, levels=c('Krandom','Krandom*N/2','ED','DTW', 'MI', 'DC'))

library('extrafont')
library(ggplot2)
p4 = ggplot(data = ggdatadens, mapping = aes(x = Time, y = AUC,color = Method,linetype = Type))+
  geom_point() + geom_line()+ggtitle("C. Time-varying topology (win=5)\n(1) without outliers") +
  theme(axis.title.x = element_text(size = 12,family = 'sans'), axis.text.x = element_text(angle = 45, hjust=1, 
                                                                                           size = 12,family = 'sans'),
        axis.title.y = element_text(size = 12,family = 'sans'), axis.text.y = element_text(size = 12,family = 'sans'),
        title=element_text(size=12,family = 'sans'),legend.text = element_text(size = 12,family = 'sans'),
        legend.position = 'bottom')+   #
  guides(color = 'none')+
  facet_grid(.~Type)+
  scale_color_manual(values=color2)
p4




