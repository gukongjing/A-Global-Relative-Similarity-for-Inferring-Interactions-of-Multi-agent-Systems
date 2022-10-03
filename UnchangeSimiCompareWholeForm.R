rm(list = ls())
library(parallel)
#Calculate the number of cores检查电脑当前可用核数
no_cores<-detectCores() - 2		#F-物理CPU核心数/T-逻辑CPU核心数 logical=F

#Initiate cluster发起集群,同时创建数个R进行并行计算
#只是创建待用的核,而不是并行运算环境
cl<-makeCluster(no_cores)

p = seq(0,0.15,0.01)
#现只需要使用并行化版本的lapply,parLapply就可以
final = parLapply(cl, p,function(p){
  
  rept = 50
  
  Nt = 40   #数量
  uav = c(1:Nt)   #无人机索引
  td = sample(1:3, 1, replace = FALSE)
  N = 300
  
  
  auc = rep(0,rept*5*2)
  dim(auc) = c(rept,5,2)
  dens = rep(0,rept)
  i=1
  for (i in 1:rept) {
    set.seed(200*i*p+1000)
    
    traj = function(k,td,Nt,N,uav,A,x0,v0){
      library(ggplot2)
      library(network)
      library(GGally)
      library(sna)
      library(RColorBrewer)
      # td =1 
      t <- 0.1
      gnum <- Nt/k
      R <- 1
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
      
      # 加速度
      a <- rep(0,Nt*N*2)
      ag <- rep(0,Nt*N*2)    # 加高斯
      dim(a) <- c(Nt,N,2)
      dim(ag) <- c(Nt,N,2)
      
      a[,1:(N-1),] = v[,2:N,]-v[,1:(N-1),]
      a[,N,] = a[,N-1,]
      i = 0
      for (i in 1:N) {
        temp = sapply(uav, function(uav){
          a[uav,i,] + rnorm(2, mean = 0, sd = sdg)
        })
        ag[,i,] = t(temp)
      }
      
      
      # 归一化
      i = 1
      for (i in 1:2) {
        xg[,,i] = (xg[,,i]-mean(xg[,,i]))/sqrt(sum((xg[,,i]-mean(xg[,,i]))^2)/(length(xg[,,i])-1))  # z-score 标准化
      }
      i = 1
      for (i in 1:2) {
        vg[,,i] = (vg[,,i]-mean(vg[,,i]))/sqrt(sum((vg[,,i]-mean(vg[,,i]))^2)/(length(vg[,,i])-1))  # z-score 标准化
      }
      i = 1
      for (i in 1:2) {
        ag[,,i] = (ag[,,i]-mean(ag[,,i]))/sqrt(sum((ag[,,i]-mean(ag[,,i]))^2)/(length(ag[,,i])-1))  # z-score 标准化
      }
      
      resu= list(xg,vg,ag)
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
    
    ### 5 mutual information
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
    
    ### 6 Distance correlation
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
    
    library(igraph)
    library(network)
    library(minet)
    
    # 生成参考网络
    ## 基础矩阵:类对角矩阵
    temp = Nt-1
    Atemp = diag(temp)
    A = matrix(0,nrow = Nt, ncol = Nt)     #类对角
    A[2:Nt,1:(Nt-1)] = Atemp
    A[upper.tri(A)] = t(A)[upper.tri(t(A))]
    LDJ = A     #类对角
    ## ER network
    ernet = erdos.renyi.game(Nt, p, type=c("gnp"), directed = FALSE, loops = FALSE)
    # plot(ernet)
    ernet = as.matrix(get.adjacency(ernet))
    gnet = ernet + LDJ
    gnet[which(gnet>1)] = 1
    diag(gnet)=0
    
    dens[i] = sum(gnet)/(Nt*(Nt-1))
    
    # G = graph.adjacency(gnet,mode="undirected")
    # plot(G)
    
    # 生成轨迹信息
    x0 <- runif(Nt*2*(1+td),min = 0,max = 10)
    v0 <- runif(Nt*2*(1+td),min = 0,max = 1)
    dim(x0) = c(Nt,(td+1),2)
    dim(v0) = c(Nt,(1+td),2)
    xvg = traj(k=2,td,Nt,N,uav,gnet,x0,v0)
    
    plot(xvg[[1]][Nt,,1],xvg[[1]][Nt,,2],type = "l") #,xlim=c(-20,20),ylim=c(-20,20)
    for (i0 in 1:(Nt-1)) {
      lines(xvg[[1]][i0,,1],xvg[[1]][i0,,2])
    }
    
    # 计算群体稳定性
    CV = sum(abs(sapply(uav,function(uav){
      sd(xvg[[2]][uav,,1])+sd(xvg[[2]][uav,,2])
    })))
    
    
    sequence = c(1:N)
    
    series1x = lapply(uav, function(uav){
      as.matrix(data.frame(coordsx1 = xvg[[1]][uav,1:N,1], coordsx2 = xvg[[1]][uav,1:N,2]))
    })
    series1v = lapply(uav, function(uav){
      as.matrix(data.frame(coordsx1 = xvg[[2]][uav,1:N,1], coordsx2 = xvg[[2]][uav,1:N,2]))
    })
    series1xv = lapply(uav, function(uav){
      as.matrix(data.frame(coordsx1 = xvg[[1]][uav,1:N,1], coordsx2 = xvg[[1]][uav,1:N,2],
                           coordsv1 = xvg[[2]][uav,1:N,1],coordsv2 = xvg[[2]][uav,1:N,2]))
    })
    # series2x = list()    # with outliers
    # series2xv = list()    # with outliers
    series = series1xv
    
    
    library(ROCR)
    
    ## 1 ED
    edsimilarity = EDsimi(series)
    # corrplot(edsimilarity)
    
    pred = prediction(as.vector(edsimilarity[upper.tri(edsimilarity)]),as.vector(gnet[upper.tri(gnet)]))
    pr.perf = performance(pred, measure = "aucpr")
    auc[i,1,1] = unlist(pr.perf@y.values)
    auc.perf = performance(pred, measure = "auc")
    auc[i,1,2] = unlist(auc.perf@y.values)
    
    #2 dtw
    dtwsimilarity = DTWsimi(series)
    
    pred = prediction(as.vector(dtwsimilarity[upper.tri(dtwsimilarity)]),as.vector(gnet[upper.tri(gnet)]))
    pr.perf = performance(pred, measure = "aucpr")
    auc[i,2,1] = unlist(pr.perf@y.values)
    auc.perf = performance(pred, measure = "auc")
    auc[i,2,2] = unlist(auc.perf@y.values)
    
    #4 Ksimi
    w=1  #historical impact
    Ksimilarity = Krandom(w,series)
    
    pred = prediction(as.vector(Ksimilarity[upper.tri(Ksimilarity)]),as.vector(gnet[upper.tri(gnet)]))
    pr.perf = performance(pred, measure = "aucpr")
    auc[i,3,1] = unlist(pr.perf@y.values)
    auc.perf = performance(pred, measure = "auc")
    auc[i,3,2] = unlist(auc.perf@y.values)
    
    #5 MI
    MIsimilarity = MIsimi(series)
    
    pred = prediction(as.vector(MIsimilarity[upper.tri(MIsimilarity)]),as.vector(gnet[upper.tri(gnet)]))
    pr.perf = performance(pred, measure = "aucpr")
    auc[i,4,1] = unlist(pr.perf@y.values)
    auc.perf = performance(pred, measure = "auc")
    auc[i,4,2] = unlist(auc.perf@y.values)
    
    #6 pcov
    pcorsimilarity = DCsimi(series)
    
    pred = prediction(as.vector(pcorsimilarity[upper.tri(pcorsimilarity)]),as.vector(gnet[upper.tri(gnet)]))
    pr.perf = performance(pred, measure = "aucpr")
    auc[i,5,1] = unlist(pr.perf@y.values)
    auc.perf = performance(pred, measure = "auc")
    auc[i,5,2] = unlist(auc.perf@y.values)
  }
  
prauc = apply(auc[,,1], 2, mean)
rocauc = apply(auc[,,2], 2, mean)
dens = mean(dens)

D = c(dens,prauc,rocauc,CV)

return(D)
})

stopCluster(cl)

D = do.call(rbind, final)
D = as.matrix(D)

# write.table(D,file = 'Dnormal0324.txt', sep = ' ', row.names = FALSE, col.names = FALSE)
D = as.matrix(read.table('Dnormal0324.txt',sep = ' ', header = FALSE))

dens = D[,1]
prauc = D[,2:6]
rocauc = D[,7:11]
CV = D[,12]

library(RColorBrewer)
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
color1 <- colorRampPalette(rev(brewer.pal(11,'Spectral')))(28)
# color2 = c("#dc745c","#efaa6e","#f2c368","#6a84b5","#3b5899")   #"#b0bcd6",
color2 = c("#dc745c","#b0bcd6","#6a84b5","#3b5899",'#041656')  #"#f2c368",

# 画图AUC-Netdens
ggdatadens = data.frame(AUC = c(as.vector(prauc),as.vector(rocauc)), 
                      Method = rep(rep(c('ED','DTW','Krandom','MI','DC'),each = length(p)),2),
                      Density = rep(rep(dens, 5),2),
                      Type = rep(c("PR",'ROC'),each = length(p)*5))
ggdatadens$Method <- factor(ggdatadens$Method, levels=c('Krandom', 'ED', 'DC','DTW', 'MI'))

library('extrafont')
library(ggplot2)
p1 = ggplot(data = ggdatadens, mapping = aes(x = Density, y = AUC,color = Method,linetype = Type))+
  geom_point() + geom_line()+ggtitle("A. Fixed topology (win=P=300)\n(1) without outliers") +
  theme(axis.title.x = element_text(size = 12,family = 'sans'), axis.text.x = element_text(angle = 45, hjust=1, 
                                                                                           size = 12,family = 'sans'),
        axis.title.y = element_text(size = 12,family = 'sans'), axis.text.y = element_text(size = 12,family = 'sans'),
        title=element_text(size=12,family = 'sans'),legend.position = 'none')+ 
  scale_color_manual(values=color2)
p1

# # 画图AUC-CV
# ggdata = data.frame(AUC = c(as.vector(prauc),as.vector(rocauc)), 
#                       Method = rep(rep(c('ED','DTW','Spca','Krandom','MI','DC'),each = length(p)),2),
#                       CV = rep(rep(CV, 6),2),
#                       Type = rep(c("PR",'ROC'),each = length(p)*6))
# library(ggplot2)
# p2 = ggplot(data = ggdata, mapping = aes(x = CV, y = AUC,color = Method,linetype = Type))+
#   geom_point() + geom_line()+
#   theme(axis.title.x = element_text(size = 12), axis.text.x = element_text(angle = 45, hjust=1, size = 12),
#         axis.title.y = element_text(size = 12), axis.text.y = element_text(size = 12),
#         title=element_text(size=15),legend.text = element_text(size = 12))+ 
#   scale_color_manual(values=color2)
# p2


# library(cowplot)
# gg1 <- ggdraw() + 
#   draw_plot(p2, 0,0,1, 1) + draw_plot(p1, 0.4,0.48,0.6,0.5)
# print(gg1)
# 
# 
# 
# #######################重新画图#########################
# D = as.matrix(read.table('D.txt',sep = ' ', header = FALSE))
# dens = D[,1]
# prauc = D[,2:6]
# rocauc = D[,7:11]
# CV = D[,12]
# 
# library(RColorBrewer)
# rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
# color1 <- colorRampPalette(rev(brewer.pal(11,'Spectral')))(28)
# color2 = c("#dc745c","#efaa6e","#f2c368","#6a84b5","#3b5899")   #"#b0bcd6",
# 
# # 画图PRAUC
# ggdatadens = data.frame(AUC = as.vector(prauc), 
#                         Method = rep(c('ED','DTW','Krandom','MI','DC'),each = length(p)),
#                         Density = rep(dens, 5))
# ggdatadens$Method <- factor(ggdatadens$Method, levels=c('Krandom', 'ED', 'DC','DTW', 'MI'))
# 
# library('extrafont')
# library(ggplot2)
# p1 = ggplot(data = ggdatadens, mapping = aes(x = Density, y = AUC,color = Method))+
#   geom_point() + geom_line()+
#   theme(axis.title.x = element_text(size = 12,family = 'sans'), axis.text.x = element_text(angle = 45, hjust=1, 
#                                                                                            size = 12,family = 'sans'),
#         axis.title.y = element_text(size = 12,family = 'sans'), axis.text.y = element_text(size = 12,family = 'sans'),
#         title=element_text(size=15,family = 'sans'),legend.position = 'none')+ 
#   scale_color_manual(values=color2)
# p1
# 
# # 画图ROCAUC
# ggdata = data.frame(AUC = as.vector(rocauc),
#                     Method = rep(c('ED','DTW','Krandom','MI','DC'),each = length(p)),
#                     Density = rep(dens, 5))
# library(ggplot2)
# p2 = ggplot(data = ggdata, mapping = aes(x = Density, y = AUC,color = Method))+
#   geom_point() + geom_line()+
#   theme(axis.title.x = element_text(size = 12), axis.text.x = element_text(angle = 45, hjust=1, size = 12),
#         axis.title.y = element_text(size = 12), axis.text.y = element_text(size = 12),
#         title=element_text(size=15),legend.text = element_text(size = 12))+
#   scale_color_manual(values=color2)
# p2
# 
# library(cowplot)
# gg1 <- ggdraw() + 
#   draw_plot(p1, 0,0,0.4, 1) + draw_plot(p2, 0.4,0,0.6,1)
# print(gg1)









