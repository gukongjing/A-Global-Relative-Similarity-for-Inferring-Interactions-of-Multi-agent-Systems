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
  
  netconstruction  = function(Nt,kind,disturbnum,k,m){
    ###构建10种不同的邻接矩阵
    
    ###输入
    #【1】Nt: 集群规模,integer，范围5-500.
    #【2】kind：指定矩阵类型，integer，范围1-10.其中1代表BCJ不对称层级矩阵,2代表CJ对称层级矩阵,3代表BLDJ不对称类对角,4代表LDJ对称类对角,
    #5代表BDC不对称矩阵,6代表DC对称矩阵,7代表BFK分块不对称矩阵,8代表FK分块对称矩阵,9代表SSJ对称上三角矩阵,10代表FSSJ分块上三角矩阵
    #【3】disturbnum: 矩阵随机连边数，integer，范围0-Nt，矩阵5 6 7 8 9 10都需要该参数。默认值设为0.1*Nt
    #【4】k：分块数目，integer，范围3-10，矩阵7 8 10需要该参数。默认值设为5
    #【5】m：分支数，integer，范围3-8，矩阵1 2需要该参数。默认值为5
    
    ###输出
    #【1】A：邻接矩阵，Nt*Nt，元素为0或1
    
    #默认设置
    # Nt = 150     
    # disturbnum = 0.1*Nt    
    # k = 5         
    # m = 5         
    
    #其他需要计算的参数
    node = c(1:Nt)     #点集
    gnum <- Nt/k   #每块数目
    n = round(log((Nt*(m-1)+1),m))  #计算层级结构矩阵的层级数
    
    ##1 层级结构
    scale = (1-m^n)/(1-m)
    A = matrix(0, nrow = max(scale,Nt), ncol = max(scale,Nt))
    for (i in 1:((1-m^(n-1))/(1-m))){
      A[i,((i-1)*m+2):(i*m+1)] = 1
    }
    if(Nt>scale){
      for (i in (scale+1):Nt) {
        temp = sample(node[-i],1)
        A[i,temp] = 1
      }
    }else{
      A = A[1:Nt,1:Nt]
    }
    CJ = A
    
    
    ##2  基础矩阵:类对角矩阵
    temp = Nt-1
    Atemp = diag(temp)
    A = matrix(0,nrow = Nt, ncol = Nt)     #类对角
    A[2:Nt,1:(Nt-1)] = Atemp
    A[upper.tri(A)] = t(A)[upper.tri(t(A))]
    LDJ = A     #类对角
    
    ##3  分块矩阵
    C = matrix(0,nrow=Nt,ncol=Nt)
    for (a in 0:(k-1)) {
      C[(a*gnum+1):((a+1)*gnum),(a*gnum+1):((a+1)*gnum)] = 1
      C[(a*gnum+1),(gnum*(a+1)+1)%%Nt] = 1
      C[(gnum*(a+1)+1)%%Nt,(a*gnum+1)] = 1
    }
    diag(C)=0
    FK = C     #分块对称矩阵
    
    
    ##4  ER network
    ernet = erdos.renyi.game(Nt, 0.01, type=c("gnp"), directed = FALSE, loops = FALSE)
    # plot(ernet)
    ernet = as.matrix(get.adjacency(ernet))
    ernet = ernet + LDJ
    
    ## ba network
    banet = sample_pa(Nt, power = 0.9, m = NULL, out.dist = NULL, out.seq = NULL,
                      out.pref = FALSE, zero.appeal = 1, directed = FALSE,
                      algorithm = c("psumtree", "psumtree-multiple", "bag"),
                      start.graph = NULL)
    # plot(banet)
    banet = as.matrix(get.adjacency(banet))
    
    # ws<-watts.strogatz.game(1,Nt,2,0.5)
    # plot(ws)
    
    #输出
    if(kind == 1){
      return(CJ)
    }else if(kind == 2){
      return(LDJ)
    }else if(kind == 3){
      return(FK)
    }else if(kind == 4){
      return(ernet)
    }else if(kind == 5){
      return(banet)
    }
    
  }
  
  
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
  
  
  ### 4 随机K   多变量
  KrandomNN = function(w,series){    #k=clustering number, w=historical impact,nd=net distance
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
  
  
  Nt <- 40   #数量
  uav = c(1:Nt)   #无人机索引
  kind = 3        # 网络类型
  k1 = 5           #子群数目
  k2 = 8
  m = 3
  disturbnum = 400
  
  N1 = 100
  N2 = 100
  N3 = 101
  N = N1+N2+N3
  td = 1

  CPreal = c(N1,N1+N2)
  
  # 生成参考网络
  g1 = netconstruction(Nt,kind,disturbnum,k1,m)
  g1[which(g1>1)] = 1
  diag(g1)=0
  G = graph.adjacency(g1,mode="undirected")
  plot(G)
  
  # 生成轨迹信息
  x0 <- runif(Nt*2*(1+td),min = 0,max = 10)
  v0 <- runif(Nt*2*(1+td),min = 0,max = 1)
  dim(x0) = c(Nt,(td+1),2)
  dim(v0) = c(Nt,(1+td),2)
  xvg1 = traj(k1,td,Nt,N1,uav,g1,x0,v0,R = 1)
  
  
  ### second stage
  # 生成参考网络
  g2 = netconstruction(Nt,kind,disturbnum,k2,m)
  g2[which(g2>1)] = 1
  diag(g2)=0
  G = graph.adjacency(g2,mode="undirected")
  plot(G)
  
  ## 轨迹
  x0 = xvg1[[1]][,((N1-td):N1),]
  v0 = xvg1[[2]][,((N1-td):N1),]
  # plot(x0[,1,1],x0[,1,2])
  # plot(x0[,2,1],x0[,2,2])
  # plot(v0[,1,1],v0[,1,2])
  
  
  xvg2 = traj(k2,td,Nt,N2+td+1,uav,g2,x0,v0,R=1)
  
  ### third stage
  # 生成参考网络
  g3 = netconstruction(Nt,kind=4,disturbnum,k2,m)
  g3[which(g3>1)] = 1
  diag(g3)=0
  G = graph.adjacency(g3,mode="undirected")
  plot(G)
  
  
  ## 轨迹
  x0 = xvg2[[1]][,((N2-td):(N2)),]
  v0 = xvg2[[2]][,((N2-td):(N2)),]
  # plot(x0[,1,1],x0[,1,2])
  # plot(x0[,2,1],x0[,2,2])
  # plot(v0[,1,1],v0[,1,2])
  
  
  xvg3 = traj(k2,td,Nt,N3+td+1,uav,g3,x0,v0,R=1)
  
  xg = rep(0,Nt*N*2)
  dim(xg) = c(Nt,N,2)
  xg[,1:N1,] = xvg1[[1]]
  xg[,(N1-td):(N1+N2),] = xvg2[[1]]
  xg[,(N1+N2-td):N,] = xvg3[[1]]
  
  
  vg = rep(0,Nt*N*2)
  dim(vg) = c(Nt,N,2)
  vg[,1:N1,] = xvg1[[2]]
  vg[,(N1-td):(N1+N2),] = xvg2[[2]]
  vg[,(N1+N2-td):N,] = xvg3[[2]]
  
  plot(xg[Nt,,1],xg[Nt,,2],type = "l") #,xlim=c(-20,20),ylim=c(-20,20)
  for (i0 in 1:(Nt-1)) {
    lines(xg[i0,,1],xg[i0,,2])
  }
  
  l = 5
  sequence = c(1:l)
  
  start = 6
  identification = c(start:(N-1))
  auc = rep(0,length(identification)*2*2)
  dim(auc) = c(length(identification),2,2)
  i = 0
  
  CV = vector()    #记录速度方差
  
  for (time in identification) {
    i = i+1
    temp = c((time-l):(time-1))
    
    CV[i] = mean(sqrt((vg[,temp,1]-mean(vg[,temp,1]))^2+(vg[,temp,2]-mean(vg[,temp,2]))^2))
    
    xstand1 = (xg[,temp,1]-mean(xg[,temp,1]))/sqrt(sum((xg[,temp,1]-mean(xg[,temp,1]))^2)/(length(xg[,temp,1])-1))
    vstand1 = (vg[,temp,1]-mean(vg[,temp,1]))/sqrt(sum((vg[,temp,1]-mean(vg[,temp,1]))^2)/(length(vg[,temp,1])-1))
    xstand2 = (xg[,temp,2]-mean(xg[,temp,2]))/sqrt(sum((xg[,temp,2]-mean(xg[,temp,2]))^2)/(length(xg[,temp,2])-1))
    vstand2 = (vg[,temp,2]-mean(vg[,temp,2]))/sqrt(sum((vg[,temp,2]-mean(vg[,temp,2]))^2)/(length(vg[,temp,2])-1))
    
    series1xv = lapply(uav, function(uav){
      as.matrix(data.frame(coordsx1 = xstand1[uav,], coordsx2 = xstand2[uav,],
                           coordsv1 = vstand1[uav,], coordsv2 = vstand2[uav,]))
    })
    series = series1xv
    
    if(time<=N1){
      gnet = g1
    }else if(time<=N1+N2){
      gnet = g2
    }else{
      gnet = g3
    }
    
    library(ROCR)
    
    #4 Ksimi
    w=1  #historical impact
    Ksimilarity = Krandom(w,series)
    pred = prediction(c(as.vector(Ksimilarity[upper.tri(Ksimilarity)]),0),
                      c(as.vector(gnet[upper.tri(gnet)]),0))
    pr.perf = performance(pred, measure = "aucpr")
    auc[i,1,1] = unlist(pr.perf@y.values)
    auc.perf = performance(pred, measure = "auc")
    auc[i,1,2] = unlist(auc.perf@y.values)
    
    Ksimilarity = KrandomNN(w,series)
    pred = prediction(c(as.vector(Ksimilarity[upper.tri(Ksimilarity)]),0),
                      c(as.vector(gnet[upper.tri(gnet)]),0))
    pr.perf = performance(pred, measure = "aucpr")
    auc[i,2,1] = unlist(pr.perf@y.values)
    auc.perf = performance(pred, measure = "auc")
    auc[i,2,2] = unlist(auc.perf@y.values)
  }
  
  prauc = auc[,,1]
  rocauc = auc[,,2]
  
  D = c(as.vector(prauc),as.vector(rocauc),CV)
  
  return(D)
})

stopCluster(cl)

D = do.call(rbind, final)

# write.table(D,file = 'switch_comparison.txt', sep = ' ', row.names = FALSE, col.names = FALSE)
D = as.matrix(read.table('switch_comparison.txt',sep = ' ', header = FALSE))

E = apply(D, 2, mean)
b = length(E)

prauc = E[1:(b/5*2)]
rocauc = E[(b/5*2+1):(b/5*4)]
VV = E[(b/5*4+1):b]
time = c(6:300)

library(RColorBrewer)
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
color1 <- colorRampPalette(rev(brewer.pal(11,'Spectral')))(28)
color2 = c("#dc745c","#f2c368","black",'#707980',"#3b5899",'#041656','black',"#b0bcd6")  #"#efaa6e",,



# 画图AUC-Netdens
ggdatadens = data.frame(Value = E, 
                        Time = rep(c(time),5),
                        Type = c(rep('PR',2*length(time)),rep('ROC',2*length(time)),
                                 rep('Velocity Variance',length(time))),
                        Method = c(rep(rep(c('Krandom','Krandom*N/2'),each = length(time)), 2),
                                   rep('Velocity Variance',length(time))))
# ,
# face = c(rep('AUC',2*length(time)),rep('Convergence',length(time)))

library('extrafont')
library(ggplot2)
gg2 = ggplot(data = ggdatadens, mapping = aes(x = Time, y = Value,color = Method))+
  geom_line()+ggtitle("B. Switching topology (win=5)\n     without outliers") +
  theme(axis.title.x = element_text(size = 12,family = 'sans'), axis.text.x = element_text(angle = 45, hjust=1, 
                                                                                           size = 12,family = 'sans'),
        axis.title.y = element_text(size = 12,family = 'sans'), axis.text.y = element_text(size = 12,family = 'sans'),
        title=element_text(size=12,family = 'sans'),legend.position = 'none')+ 
  facet_grid(Type~.,scales = "free")+
  geom_vline(aes(xintercept=100),alpha=0.4)+
  geom_vline(aes(xintercept=200),alpha=0.4)+
  scale_color_manual(values=color2)
gg2








