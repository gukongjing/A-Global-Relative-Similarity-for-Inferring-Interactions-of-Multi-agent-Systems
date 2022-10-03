rm(list = ls())
library(parallel)
#Calculate the number of cores检查电脑当前可用核数
no_cores<-detectCores() -2		#F-物理CPU核心数/T-逻辑CPU核心数 logical=F

#Initiate cluster发起集群,同时创建数个R进行并行计算
#只是创建待用的核,而不是并行运算环境
cl<-makeCluster(no_cores)

p = c(0,0.1,0.2,0.4,0.6,0.8,1)
#现只需要使用并行化版本的lapply,parLapply就可以
final = parLapply(cl, p,function(p){
  library(forcats)
  
  rept = 2
  traj = function(k,td,Nt,N,uav,A,x0,v0,inside){
    library(ggplot2)
    library(network)
    library(GGally)
    library(sna)
    library(RColorBrewer)
    # td =1 
    t <- 0.1
    # gnum <- Nt/k
    R <- 1
    td = td     # time delay
    sdg = 0.01    # 高斯噪声的方差   
    out = p    # p比例错误
    
    
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
    xg1 <- rep(0,Nt*N*2)    # 加高斯和离散异常值
    vg1 <- rep(0,Nt*N*2)    # 加高斯和离散异常值
    xg2 <- rep(0,Nt*N*2)    # 加高斯和成片异常值
    vg2 <- rep(0,Nt*N*2)    # 加高斯和成片异常值
    
    dim(x) <- c(Nt,N,2)
    dim(v) <- c(Nt,N,2)
    dim(xg) <- c(Nt,N,2)    
    dim(vg) <- c(Nt,N,2)    
    dim(xg1) <- c(Nt,N,2)    
    dim(vg1) <- c(Nt,N,2)    
    dim(xg2) <- c(Nt,N,2)    
    dim(vg2) <- c(Nt,N,2)    
    
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
    
    
    plot(xg[1,25:N,1],xg[1,25:N,2],xlab = "x-coordinate",ylab = "y-coordinate",type = 'l',
         cex.lab=2,cex.axis=2)   # cex.main=1.5,
    for (z in 2:Nt) {
      lines(xg[z,25:N,1],xg[z,25:N,2])
    }
    
    xg1 = xg
    vg1 = vg
    
    # 添加高斯噪声和离散异常值
    outp = sample(c(1:N),floor(out*N))
    outnum = floor(inside*Nt*2*2) 
    pool = c(as.vector(xg),as.vector(vg))
    
    if(p>0 && inside>0){
      for (oo in outp) {
        temp = rep(0,Nt*2*2)
        dim(temp) = c(Nt,4)
        temp[,1:2] = xg[,oo,]
        temp[,3:4] = vg[,oo,]
        temp1 = as.vector(temp)
        ind = sample(length(temp1),outnum)
        temp1[ind] = pool[sample(length(pool),outnum)]
        dim(temp1) = c(Nt,4)
        xg1[,oo,] = temp1[,1:2]
        vg1[,oo,] = temp1[,3:4]
      }
    }
    
    
    sum(xg!=xg1)+ sum(vg!=vg1)
    
    i=0
    ggxg = xg[1,25:(N),]
    ggxg1 = xg1[1,25:(N),]
    for (i in 2:Nt) {
      ggxg = rbind(ggxg,xg[i,25:(N),])
      ggxg1 = rbind(ggxg1,xg1[i,25:(N),])
    }
    
    outliers1 = which(ggxg1[,1]!=ggxg[,1])
    outliers2 = which(ggxg1[,2]!=ggxg[,2])
    outliers = union(outliers1,outliers2)
    
    ggseg = data.frame(
      traj.x1 = ggxg[outliers,1],   
      traj.x2 = ggxg[outliers,2]
    )
    point = factor(rep(c(1:Nt),each = N-24))
    # 
    # # ggplot traj
    ggdata <- data.frame(
      poin = point[-outliers],
      traj.x1=ggxg1[-outliers,1],
      traj.x2=ggxg1[-outliers,2]
    )
    
    display.brewer.pal(9,"Greys")
    grey = brewer.pal(9,"Greys")
    mycolors<- grey[c(5,6,7,8,9)]
    p1 = ggplot() +  geom_path(data = ggdata, mapping = aes(x = traj.x1, y=traj.x2,group = poin),size = 0.5) +
      geom_point(data = ggseg,aes(x = traj.x1, y=traj.x2),color = '#FFC83F',alpha = 1,size = 0.5) +
      theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
            axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(),
            title=element_text(size=12),legend.text = element_text(size = 15),legend.position = c(0.85,0.25))+
      # scale_colour_manual(values=c('black','grey75','dimgray'))+
      scale_color_manual(values = mycolors)+
      ggtitle("A. Random outliers\n    (p1=0.1, p2=0.1)") +
      theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank()) +  # 去掉背景 
      guides(color='none')
    
    # 归一化
    i = 1
    for (i in 1:2) {
      xg1[,,i] = (xg1[,,i]-mean(xg1[,,i]))/sqrt(sum((xg1[,,i]-mean(xg1[,,i]))^2)/(length(xg1[,,i])-1))  
    }
    i = 1
    for (i in 1:2) {
      vg1[,,i] = (vg1[,,i]-mean(vg1[,,i]))/sqrt(sum((vg1[,,i]-mean(vg1[,,i]))^2)/(length(vg1[,,i])-1)) 
    }
    
    # 添加高斯噪声和成片异常值
    ###确定成片的位置
    outlength = floor(out*N)
    outstart = sample(c(1:(N-outlength+1)),1)
    
    xg2 = xg
    vg2 = vg
    
    if(p>0 && inside>0){
      for (oo in outstart:(outstart+outlength-1)) {
        temp = rep(0,Nt*2*2)
        dim(temp) = c(Nt,4)
        temp[,1:2] = xg[,oo,]
        temp[,3:4] = vg[,oo,]
        temp1 = as.vector(temp)
        ss = sample(length(temp1)-outnum+1,1)  
        ind = c(ss:(ss+outnum-1))
        temp1[ind] = pool[sample(length(pool),outnum)]
        dim(temp1) = c(Nt,4)
        xg2[,oo,] = temp1[,1:2]
        vg2[,oo,] = temp1[,3:4]
      }
    }
    
    
    sum(xg!=xg2)+ sum(vg!=vg2)
    
    
    i=0
    ggxg = xg[1,25:(N),]
    ggxg1 = xg2[1,25:(N),]
    for (i in 2:Nt) {
      ggxg = rbind(ggxg,xg[i,25:(N),])
      ggxg1 = rbind(ggxg1,xg2[i,25:(N),])
    }
    
    outliers1 = which(ggxg1[,1]!=ggxg[,1])
    outliers2 = which(ggxg1[,2]!=ggxg[,2])
    outliers = union(outliers1,outliers2)
    
    ggseg = data.frame(
      traj.x1 = ggxg[outliers,1],   
      traj.x2 = ggxg[outliers,2]
    )
    point = factor(rep(c(1:Nt),each = N-24))
    # 
    # # ggplot traj
    ggdata <- data.frame(
      poin = point[-outliers],
      traj.x1=ggxg1[-outliers,1],
      traj.x2=ggxg1[-outliers,2]
    )
    
    display.brewer.pal(9,"Greys")
    grey = brewer.pal(9,"Greys")
    mycolors<- grey[c(5,6,7,8,9)]
    p2 = ggplot() +  geom_path(data = ggdata, mapping = aes(x = traj.x1, y=traj.x2,group = poin),size = 0.5) +
      geom_point(data = ggseg,aes(x = traj.x1, y=traj.x2),color = '#FFC83F',alpha = 1,size = 0.5) +
      theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
            axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(),
            title=element_text(size=12),legend.text = element_text(size = 15),legend.position = c(0.85,0.25))+
      # scale_colour_manual(values=c('black','grey75','dimgray'))+
      scale_color_manual(values = mycolors)+
      ggtitle("B. Concentrated outliers\n    (p1=0.1, p2=0.1)") +
      theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank()) +  # 去掉背景 
      guides(color='none')
    
    
    # 归一化
    i = 1
    for (i in 1:2) {
      xg2[,,i] = (xg2[,,i]-mean(xg2[,,i]))/sqrt(sum((xg2[,,i]-mean(xg2[,,i]))^2)/(length(xg2[,,i])-1))  
    }
    i = 1
    for (i in 1:2) {
      vg2[,,i] = (vg2[,,i]-mean(vg2[,,i]))/sqrt(sum((vg2[,,i]-mean(vg2[,,i]))^2)/(length(vg2[,,i])-1)) 
    }
    
    plot(xg2[1,25:N,1],xg2[1,25:N,2],xlim=c(0,3),ylim=c(0,3),xlab = "x-coordinate",ylab = "y-coordinate",type = 'l',
         cex.lab=2,cex.axis=2)   # cex.main=1.5,
    for (z in 2:Nt) {
      lines(xg2[z,25:N,1],xg2[z,25:N,2])
    }
    
    resu= list(xg,vg,xg1,vg1,xg2,vg2)
    return(resu)
  }
  
  
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
  
  
  Nt = 40   #数量
  uav = c(1:Nt)   #无人机索引
  td = 1
  N = 301
  insideseq = seq(0,1,0.2)
  
  auc = rep(0,rept*4*length(insideseq))
  dim(auc) = c(rept,4,length(insideseq))
  i=1
  for (i in 1:rept) {
    
    library(igraph)
    library(network)
    
    # 生成参考网络
    ## 基础矩阵:类对角矩阵
    temp = Nt-1
    Atemp = diag(temp)
    A = matrix(0,nrow = Nt, ncol = Nt)     #类对角
    A[2:Nt,1:(Nt-1)] = Atemp
    A[upper.tri(A)] = t(A)[upper.tri(t(A))]
    LDJ = A     #类对角
    ## ER network
    ernet = erdos.renyi.game(Nt, 0.02, type=c("gnp"), directed = FALSE, loops = FALSE)
    # plot(ernet)
    ernet = as.matrix(get.adjacency(ernet))
    gnet = ernet + LDJ
    gnet[which(gnet>1)] = 1
    diag(gnet)=0
    
    
    # G = graph.adjacency(gnet,mode="undirected")
    # plot(G)
    
    # 生成轨迹信息
    x0 <- runif(Nt*2*(1+td),min = 0,max = 8)
    v0 <- runif(Nt*2*(1+td),min = 0,max = 1)
    dim(x0) = c(Nt,(td+1),2)
    dim(v0) = c(Nt,(1+td),2)
    iii = 1
    
    for (iii in 1:length(insideseq)) {
      inside = insideseq[iii]
      xvg = traj(k=2,td,Nt,N,uav,gnet,x0,v0,inside)
      
      sequence = c(1:N)
      
      series2 = lapply(uav, function(uav){
        as.matrix(data.frame(coordsx1 = xvg[[3]][uav,1:N,1], coordsx2 = xvg[[3]][uav,1:N,2],
                             coordsv1 = xvg[[4]][uav,1:N,1],coordsv2 = xvg[[4]][uav,1:N,2]))
      })
      series3 = lapply(uav, function(uav){
        as.matrix(data.frame(coordsx1 = xvg[[5]][uav,1:N,1], coordsx2 = xvg[[5]][uav,1:N,2],
                             coordsv1 = xvg[[6]][uav,1:N,1],coordsv2 = xvg[[6]][uav,1:N,2]))
      })
      
      library(ROCR)
      #4 Ksimi
      w=1  #historical impact
      Ksimilarity = Krandom(w,series2)
      pred = prediction(as.vector(Ksimilarity[upper.tri(Ksimilarity)]),as.vector(gnet[upper.tri(gnet)]))
      pr.perf = performance(pred, measure = "aucpr")
      auc[i,1,iii] = unlist(pr.perf@y.values)
      auc.perf = performance(pred, measure = "auc")
      auc[i,2,iii] = unlist(auc.perf@y.values)
      
      Ksimilarity = Krandom(w,series3)
      pred = prediction(as.vector(Ksimilarity[upper.tri(Ksimilarity)]),as.vector(gnet[upper.tri(gnet)]))
      pr.perf = performance(pred, measure = "aucpr")
      auc[i,3,iii] = unlist(pr.perf@y.values)
      auc.perf = performance(pred, measure = "auc")
      auc[i,4,iii] = unlist(auc.perf@y.values)
      
    }
    
  }
  
  prauc1 = apply(auc[,1,], 2, mean)
  prauc2 = apply(auc[,3,], 2, mean)
  rocauc1 = apply(auc[,2,], 2, mean)
  rocauc2 = apply(auc[,4,], 2, mean)
  
  D = c(prauc1,prauc2,rocauc1,rocauc2)
  
  return(D)
})

stopCluster(cl)

D = do.call(rbind, final)
# write.table(D,file = 'outlierstest7.txt', sep = ' ', row.names = FALSE, col.names = FALSE)
D = as.matrix(read.table('outlierstest7.txt',sep = ' ', header = FALSE))   #outlierstest2

prauc1 = D[,1:6]
prauc2 = D[,7:12]
rocauc1 = D[,13:18]
rocauc2 = D[,19:24]

library(RColorBrewer)
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
color1 <- colorRampPalette(rev(brewer.pal(11,'Spectral')))(28)
color2 = c("#dc745c","#efaa6e","#f2c368","#b0bcd6","#6a84b5","#3b5899",'#041656')

# 画图AUC-Netdens
ggdatadens = data.frame(AUC = c(as.vector(prauc1),as.vector(prauc2),as.vector(rocauc1),as.vector(rocauc2)), 
                        Method = rep(rep(rep(c('Random','Concentrated'),
                                         each = length(p)),each = 6),2),
                        p1 = rep(rep(p, 6),4),
                        p2 = as.factor(rep(rep(seq(0,1,0.2),each = length(p)),4)),
                        Type = rep(c("PR",'ROC'),each = length(p)*6*2))
ggdatadens$Method <- factor(ggdatadens$Method, levels=c('Random','Concentrated'))

library('extrafont')
library(ggplot2)
p3 = ggplot(data = ggdatadens, mapping = aes(x = p1, y = AUC,color = p2,linetype = Type))+
  geom_point() + geom_line()+ggtitle("C. Proportion of snapshots - AUC") +
  theme(axis.title.x = element_text(size = 12,family = 'sans'), axis.text.x = element_text(angle = 45, hjust=1, 
                                                                                           size = 12,family = 'sans'),
        axis.title.y = element_text(size = 12,family = 'sans'), axis.text.y = element_text(size = 12,family = 'sans'),
        title=element_text(size=12,family = 'sans'),legend.text = element_text(size = 12,family = 'sans'))+ 
  facet_grid(.~Method)+
  scale_color_manual(values=color2)
p3


#############################################
color2 = c("#dc745c","#efaa6e","#f2c368",'#d8e1e8',"#b0bcd6","#6a84b5","#3b5899",'#041656',
           "#6a84b6","#3b5810",'#041657')
ggdatadens = data.frame(AUC = c(as.vector(prauc1),as.vector(prauc2),as.vector(rocauc1),as.vector(rocauc2)), 
                        Method = rep(rep(rep(c('Random','Concentrated'),
                                             each = length(p)),each = 6),2),
                        p1 = as.factor(rep(rep(p, 6),4)),
                        p2 = rep(rep(seq(0,1,0.2),each = length(p)),4),
                        Type = rep(c("PR",'ROC'),each = length(p)*6*2))
ggdatadens$Method <- factor(ggdatadens$Method, levels=c('Random','Concentrated'))

library('extrafont')
library(ggplot2)
p4 = ggplot(data = ggdatadens, mapping = aes(x = p2, y = AUC,color = p1,linetype = Type))+
  geom_point() + geom_line()+ggtitle("D. Proportion of outliers on a snapshot - AUC") +
  theme(axis.title.x = element_text(size = 12,family = 'sans'), axis.text.x = element_text(angle = 45, hjust=1, 
                                                                                           size = 12,family = 'sans'),
        axis.title.y = element_text(size = 12,family = 'sans'), axis.text.y = element_text(size = 12,family = 'sans'),
        title=element_text(size=12,family = 'sans'),legend.text = element_text(size = 12,family = 'sans'))+ 
  facet_grid(.~Method)+
  scale_color_manual(values=color2)
p4

library(cowplot)
gg1 <- ggdraw() + 
  draw_plot(p1, 0,0.5,0.22, 0.5) + draw_plot(p2, 0,0,0.22,0.5)+ draw_plot(p3, 0.22,0,0.39,1)+
  draw_plot(p4, 0.61,0,0.39,1)
print(gg1)

