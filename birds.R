rm(list = ls(all = TRUE))

####similarity measurement
library(dtwclust)
library(dtw)
library(forcats)
library(igraph)
library(SimilarityMeasures)
library(infotheo)
library(ggm)
library(energy)
library(segmented)
library(aricode)

Threshold = function(C){    #k=clustering number, w=historical impact,nd=net distance
  library(tidyverse) # 加载tidyverse包
  library(mclust) # 加载mclust包
  C0 = as.vector(C/max(C))
  mod5 <- Mclust(C0,G=2)
  # plot(mod5, what = "classification")
  summary(mod5)
  
  C0 = rbind(C0,mod5[["classification"]])
  C1 = C0[,order(C0[1,])]
  C1 = C1[,which(C1[1,]<=mean(C1[1,which(C1[2,]==2)]))]   #去掉比第二个均值大的数
  C1 = C1[,which(C1[1,]>=mean(C1[1,which(C1[2,]==1)]))]   #去掉比第一个均值小的数
  thresh = (max(C1[1,which(C1[2,]==1)])+min(C1[1,which(C1[2,]==2)]))/2    #属于类别1的最大值与属于类别2的最小值取平均
  
  d = C
  d[which(d<thresh)] = 0
  
  # corrplot(d)
  return(list(d,thresh))
}

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


##PR and ROC curve
Tcurve = function(simimatrix,gnet){
  threshseq = seq(0,1,0.02)
  curve1 = matrix(0, nrow = length(threshseq), ncol = 2)    #PR curve
  curve2 = matrix(0, nrow = length(threshseq), ncol = 2)    #ROC curve
  curveARI = matrix(0, nrow = length(threshseq), ncol = 2)    #ARI curve
  curveARI[,1] = threshseq
  i = 0
  for (thresh in threshseq) {
    i = i+1
    kk = simimatrix
    kk[which(kk<thresh)] = 0
    # corrplot(kk)
    kk[which(kk>0)] = 1
    diag(kk) = 0
    
    #计算精度
    TP = sum(kk*gnet)/2
    TN = (sum((1-kk)*(1-gnet))-Nt)/2
    FP = sum(kk*(1-gnet))/2
    FN = sum((1-kk)*gnet)/2
    precision=TP/(TP+FP)
    recall = TP/(TP+FN)
    curve1[i,] = c(recall,precision) 
    
    FPR = FP/(FP+TN)
    TPR = TP/(TP+FN)
    curve2[i,] = c(FPR,TPR)
    
    curveARI[i,2] = ARIvalue(kk,stand)
  }
  ind11 = order(curve1[,2])
  curve1 = curve1[ind11,]
  ind12 = order(curve1[,1])
  curve1 = curve1[ind12,]
  
  ind21 = order(curve2[,2])
  curve2 = curve2[ind21,]
  ind22 = order(curve2[,1])
  curve2 = curve2[ind22,]
  # plot(curve2)
  
  curve11 = c(curve1[-1,1],1)
  curve12 = c(0,curve1[-1,1])
  PRauc = curve1[,2]*(curve11-curve12)
  PRauc = sum(PRauc)
  
  curve21 = c(curve2[-1,1],1)
  curve22 = c(0,curve2[-1,1])
  ROCauc = curve2[,2]*(curve21-curve22)
  ROCauc = sum(ROCauc)
  
  return(list(curve1,curveARI,curve2))
}


###############暂时不用截断
ARIvalue = function(mat,gnet){
  a = mat[upper.tri(mat)]
  b = gnet[upper.tri(gnet)]
  return(ARI(a,b))
}


library(igraph)
library(forcats)
library(ROCR)
library(aricode)
library(energy)

Nt = 9
uav = seq(1:Nt)

# standard answer
hierarchy = c(2,9,6,7,5,1,4,8,3)           # maximum means highest hierarchy
names(hierarchy) = c(1:Nt)
barplot(hierarchy,ylab = "ranking",main = "hierarchy",sub = "maximum means highest hierarchy")

delay = matrix(0, nrow = Nt, ncol = Nt)
delay[1,] = c(1,0.6,1,0.2,1,1,1,0.2,1)
delay[2,] = c(1,1,1,1,1,1,1,1,1)
delay[3,] = c(1,0.6,1,0.2,1,1,1,0.2,1)
delay[4,] = c(1,0.4,1,1,1,1,1,0.2,1)
delay[5,] = c(1,0.4,1,0.2,1,1,1,0.2,1)
delay[6,] = c(1,0.8,0.2,0.4,0.2,1,0.2,0.4,0.2)
delay[7,] = c(1,0.6,0.2,0.2,1,1,1,0.4,1)
delay[8,] = c(1,0.4,1,1,1,1,1,1,1)
delay[9,] = c(1,0.6,1,0.2,1,1,1,0.4,1)

gnet = 1-delay
gnet.undirected = gnet+t(gnet)
gtemp = gnet
gtemp[which(gtemp != 0)] =1
g1 = graph_from_adjacency_matrix(gtemp,mode = 'directed')
plot(g1)
g2 = graph_from_adjacency_matrix(gtemp,mode = 'undirected')
plot(g2)
gplain = as.matrix(get.adjacency(g2,type = 'both'))      #标准答案

stand = gnet.undirected
stand[which(stand>0)] = 1

# data features: 5 times a second—— 200ms a time
# import data
filame = c("BEAB_A21_07","BEAB_A30_07","BEAB_A31_07","BEAB_A33_07",
             "BEAB_A40_07","BEAB_A46_07","BEAB_A49_07","BEAB_A63_07","BEAB_K03_07")
          
all = list()
all = lapply(uav, function(uav){
  site = paste("data\\birds\\Data_Dryad\\",filame[uav],".csv", sep = "")
  file = read.csv(site, header = T)
  data.frame(time = file$UTC.TIME,coords1 = file$LATITUDE, coords2 = file$LONGITUDE, 
             v1 = file$SPEED, v2 = file$HEADING/180)
})

######standardization
N = 11922
lat = sapply(uav,function(uav){
  all[[uav]][1:N,2]
})
lon = sapply(uav,function(uav){
  all[[uav]][1:N,3]
})
velocity = sapply(uav,function(uav){
  all[[uav]][1:N,4]
})
direct = sapply(uav,function(uav){
  all[[uav]][1:N,5]
})

lat = (lat-mean(lat))/sqrt(sum((lat-mean(lat))^2)/(length(lat)-1))  # z-score 标准化
lon = (lon-mean(lon))/sqrt(sum((lon-mean(lon))^2)/(length(lon)-1))  # z-score 标准化
velocity = (velocity-mean(velocity))/sqrt(sum((velocity-mean(velocity))^2)/(length(velocity)-1))  # z-score 标准化
direct = (direct-mean(direct))/sqrt(sum((direct-mean(direct))^2)/(length(direct)-1))  # z-score 标准化

sequence = c(1:N)

series1x = lapply(uav, function(uav){
  as.matrix(data.frame(coordsx1 = lat[1:N,uav], coordsx2 = lon[1:N,uav]))
})
series1v = lapply(uav, function(uav){
  as.matrix(data.frame(coordsx1 = velocity[1:N,uav], coordsx2 = direct[1:N,uav]))
})
series1xv = lapply(uav, function(uav){
  as.matrix(data.frame(coordsx1 = lat[1:N,uav], coordsx2 = lon[1:N,uav],
                       coordsv1 = velocity[1:N,uav],coordsv2 = direct[1:N,uav]))
})

series = series1xv




# # ggplot traj
ggdata <- data.frame(
  poin = factor(rep(c(1:Nt),each = dim(lat)[1])),
  traj.x1=as.vector(lat),
  traj.x2=as.vector(lon)
)

library(RColorBrewer)
display.brewer.pal(9,"Greys")
grey = brewer.pal(9,"Greys")
mycolors<- grey[c(1:9)]
p2 = ggplot() +  geom_path(data = ggdata, mapping = aes(x = traj.x1, y=traj.x2,color = poin, group = poin),size = 0.5) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        title=element_text(size=12),legend.text = element_text(size = 12),legend.position = c(0.85,0.25))+
  scale_color_manual(values = mycolors)+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank()) +  # 去掉背景
  guides(color='none')

threshARI = vector()
bestARI = vector()

## 1 ED
edsimilarity = EDsimi(series)
pred.ed = prediction(as.vector(edsimilarity[upper.tri(edsimilarity)]),as.vector(stand[upper.tri(stand)]))
pr.perf.ed = performance(pred.ed, measure = "aucpr")
prauc.ed = unlist(pr.perf.ed@y.values)
prcurve.ed <- performance(pred.ed, measure = "prec", x.measure = "rec") ###prcurve@x.values[[1]]为recall值

roc.perf.ed = performance(pred.ed, measure = "auc")
rocauc.ed = unlist(roc.perf.ed@y.values)
roccurve.ed <- performance(pred.ed, measure = "tpr", x.measure = "fpr")  ###roccurve@x.values[[1]]为recall值


#2 dtw
dtwsimilarity = DTWsimi(series)
pred.dtw = prediction(as.vector(dtwsimilarity[upper.tri(dtwsimilarity)]),as.vector(stand[upper.tri(stand)]))
pr.perf.dtw = performance(pred.dtw, measure = "aucpr")
prauc.dtw = unlist(pr.perf.dtw@y.values)
prcurve.dtw <- performance(pred.dtw, measure = "prec", x.measure = "rec") ###prcurve@x.values[[1]]为recall值

roc.perf.dtw = performance(pred.dtw, measure = "auc")
rocauc.dtw = unlist(roc.perf.dtw@y.values)
roccurve.dtw <- performance(pred.dtw, measure = "tpr", x.measure = "fpr")  ###roccurve@x.values[[1]]为recall值

#4 Ksimi
w=1  #historical impact
Ksimilarity = Krandom(w,series)
KThresh = Threshold(Ksimilarity)
KThreshmat = KThresh[[1]]
threshK = KThresh[[2]]
kk = KThreshmat
kk[which(kk>0)] = 1
diag(kk) = 0
library(corrplot)
corrplot(kk)

threshARI[4] = ARIvalue(kk,stand)
aa = Tcurve(Ksimilarity,gnet)
thresh4 = aa[[2]][which.max(aa[[2]][,2]),1]   #bestARI对应的threshold
bestARI = max(aa[[2]][,2])

# #####假设分为三个正态分布
# set.seed(100000000)
# mod5 <- Mclust(as.vector(Ksimilarity),G=3)
# plot(mod5, what = "density")
# C0 = rbind(as.vector(Ksimilarity),mod5[["classification"]])
# C1 = C0[,order(C0[1,])]
# C1 = C1[,which(C1[1,]<=mean(C1[1,which(C1[2,]==2)]))]   #去掉比第二个均值大的数
# C1 = C1[,which(C1[1,]>=mean(C1[1,which(C1[2,]==1)]))]   #去掉比第一个均值小的数
# thresh = (max(C1[1,which(C1[2,]==1)])+min(C1[1,which(C1[2,]==2)]))/2    #属于类别1的最大值与属于类别2的最小值取平均
# 
# d = C
# d[which(d<thresh)] = 0

pred.K = prediction(as.vector(Ksimilarity[upper.tri(Ksimilarity)]),as.vector(stand[upper.tri(stand)]))
pr.perf.K = performance(pred.K, measure = "aucpr")
prauc.K = unlist(pr.perf.K@y.values)
prcurve.K <- performance(pred.K, measure = "prec", x.measure = "rec") ###prcurve@x.values[[1]]为recall值

roc.perf.K = performance(pred.K, measure = "auc")
rocauc.K = unlist(roc.perf.K@y.values)
roccurve.K <- performance(pred.K, measure = "tpr", x.measure = "fpr")  ###roccurve@x.values[[1]]为recall值

# plot(prcurve.K)



#5 MI
MIsimilarity = MIsimi(series)
pred.MI = prediction(as.vector(MIsimilarity[upper.tri(MIsimilarity)]),as.vector(stand[upper.tri(stand)]))
pr.perf.MI = performance(pred.MI, measure = "aucpr")
prauc.MI = unlist(pr.perf.MI@y.values)
prcurve.MI <- performance(pred.MI, measure = "prec", x.measure = "rec") ###prcurve@x.values[[1]]为recall值

roc.perf.MI = performance(pred.MI, measure = "auc")
rocauc.MI = unlist(roc.perf.MI@y.values)
roccurve.MI <- performance(pred.MI, measure = "tpr", x.measure = "fpr")  ###roccurve@x.values[[1]]为recall值

#6 pcov
pcorsimilarity = DCsimi(series)
pred.pcor = prediction(as.vector(pcorsimilarity[upper.tri(pcorsimilarity)]),as.vector(stand[upper.tri(stand)]))
pr.perf.pcor = performance(pred.pcor, measure = "aucpr")
prauc.pcor = unlist(pr.perf.pcor@y.values)
prcurve.pcor <- performance(pred.pcor, measure = "prec", x.measure = "rec") ###prcurve@x.values[[1]]为recall值

roc.perf.pcor = performance(pred.pcor, measure = "auc")
rocauc.pcor = unlist(roc.perf.pcor@y.values)
roccurve.pcor <- performance(pred.pcor, measure = "tpr", x.measure = "fpr")  ###roccurve@x.values[[1]]为recall值


# 画图
library(RColorBrewer)
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
colormap <- colorRampPalette(rev(brewer.pal(11,'Spectral')))(28)
colormap = c("#dc745c","#efaa6e","#b0bcd6","#6a84b5","#3b5899")   #'#d8e1e8',"#f2c368",'#041656'
ggdata = data.frame(AUC = round(c(prauc.ed,prauc.dtw,prauc.K,prauc.MI,prauc.pcor),2), 
                    Method = c('ED','DTW','Krandom','MI','DC'))
ggdata$Method = factor(ggdata$Method,levels = c('Krandom','MI','ED','DC','DTW'))
# save(ggdata,file="data/ggdata.Rdata")
p1 = ggplot(data = ggdata, mapping = aes(x = Method, y = AUC))+   #,label = round(AUC,3)
  theme(axis.title.x = element_text(size = 12), axis.text.x = element_text(angle = 35, hjust=1, size = 12),
        axis.title.y = element_text(size = 12), axis.text.y = element_text(size = 12),
        title=element_text(size=15),legend.text = element_text(size = 12))+ 
  # geom_text(aes(y = AUC), position = position_dodge(0.9), vjust = -0.5)+
  geom_bar(stat="identity",aes(fill = Method))+guides(fill='none')+ 
  coord_cartesian(ylim = c(0.75, 0.95))+scale_y_continuous(breaks=seq(0.8, 0.9, 0.1))+
  scale_fill_manual(values=colormap)

ggdatapr = data.frame(Recall = c(prcurve.K@x.values[[1]],prcurve.ed@x.values[[1]][-1],
                                 prcurve.MI@x.values[[1]][-1], prcurve.pcor@x.values[[1]][-1],
                                 prcurve.dtw@x.values[[1]][-1]),
                      Precision = c(prcurve.K@y.values[[1]],prcurve.ed@y.values[[1]][-1],
                                    prcurve.MI@y.values[[1]][-1],prcurve.pcor@y.values[[1]][-1],
                                    prcurve.dtw@y.values[[1]][-1]),
                      Method = rep(c(paste('Krandom (AUC:',round(prauc.K,2),')'),paste('ED (AUC:',round(prauc.ed,2),')'),
                                     paste('MI (AUC:',round(prauc.MI,2),')'),
                                     paste('DC (AUC:',round(prauc.pcor,2),')'),
                                     paste('DTW (AUC:',round(prauc.dtw,2),')')), each = 36))
ggdatapr$Method = factor(ggdatapr$Method,levels = c(paste('Krandom (AUC:',round(prauc.K,2),')'),
                                                    paste('MI (AUC:',round(prauc.MI,2),')'),
                                                    paste('ED (AUC:',round(prauc.ed,2),')'),
                                                    paste('DC (AUC:',round(prauc.pcor,2),')'),
                                                    paste('DTW (AUC:',round(prauc.dtw,2),')')))
# save(ggdatapr,file="data/ggdatapr.Rdata")
p2 = ggplot(data = ggdatapr, mapping = aes(x = Recall, y = Precision,colour = Method))+
  theme(axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 12), axis.text.y = element_text(size = 12),
        title=element_text(size=12),legend.text = element_text(size = 12),
        plot.title = element_text(size = 12))+ggtitle("PR Curve") +
  # geom_text(aes(y = AUC), position = position_dodge(0.9), vjust = -0.5)+
  guides(color='none')+coord_cartesian(ylim = c(0.3, 1))+
  scale_color_manual(values=colormap)+
  theme(legend.position = c(0.7,0.7))+
  geom_line(size = 1)
  # geom_smooth(aes(x = Recall, y = Precision,colour = Method),size = 1,se = FALSE)

library(cowplot)
gg1 <- ggdraw() + 
  draw_plot(p2, 0,0,1, 1) + draw_plot(p1, 0.37,0.13,0.6,0.4)
print(gg1)



ggdata2 = data.frame(AUC = round(c(rocauc.ed,rocauc.dtw,rocauc.K,rocauc.MI,rocauc.pcor),2),
                    Method = c('ED','DTW','Krandom','MI','DC'))
ggdata2$Method = factor(ggdata2$Method,levels = c('Krandom','MI','ED','DC','DTW'))
# save(ggdata2,file="data/ggdata2.Rdata")
ggdataroc = data.frame(FPR = c(roccurve.K@x.values[[1]],roccurve.ed@x.values[[1]][-1],
                               roccurve.MI@x.values[[1]][-1],roccurve.pcor@x.values[[1]][-1],
                               roccurve.dtw@x.values[[1]][-1]),
                       TPR = c(roccurve.K@y.values[[1]],roccurve.ed@y.values[[1]][-1],
                               roccurve.MI@y.values[[1]][-1],roccurve.pcor@y.values[[1]][-1],
                               roccurve.dtw@y.values[[1]][-1]),
                       Method = rep(c(paste('Krandom (AUC:',round(rocauc.K,2),')'),
                                      paste('ED (AUC:',round(rocauc.ed,2),')'),
                                      paste('MI (AUC:',round(rocauc.MI,2),')'),
                                      paste('DC (AUC:',round(rocauc.pcor,2),')'),
                                      paste('DTW (AUC:',round(rocauc.dtw,2),')')),
                                    each = 36))
ggdataroc$Method = factor(ggdataroc$Method,levels = c(paste('Krandom (AUC:',round(rocauc.K,2),')'),
                                                      paste('MI (AUC:',round(rocauc.MI,2),')'),
                                                      paste('ED (AUC:',round(rocauc.ed,2),')'),
                                                      paste('DC (AUC:',round(rocauc.pcor,2),')'),
                                                      paste('DTW (AUC:',round(rocauc.dtw,2),')')))
# save(ggdataroc,file="data/ggdataroc.Rdata")
p3 = ggplot(data = ggdata2, mapping = aes(x = Method, y = AUC,fill = Method))+  #,label = round(AUC,3)
  theme(axis.title.x = element_text(size = 12), axis.text.x = element_text(angle = 35, hjust=1,size = 12),
        axis.title.y = element_text(size = 12), axis.text.y = element_text(size = 12),
        title=element_text(size=10),legend.text = element_text(size = 10),
        plot.title = element_text(size = 12))+ 
  # geom_text(aes(y = AUC), position = position_dodge(0.9), vjust = -0.5)+
  geom_bar(stat="identity")+guides(fill='none')+ 
  coord_cartesian(ylim = c(0.55, 0.8)) + scale_y_continuous(breaks=seq(0.6, 0.7, 0.1))+
  scale_fill_manual(values=colormap)

p4 = ggplot(data = ggdataroc, mapping = aes(x = FPR, y = TPR,colour = Method))+
  theme(axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 12), axis.text.y = element_text(size = 12),
        title=element_text(size=10),legend.text = element_text(size = 10))+ ggtitle("ROC Curve") +
  # geom_text(aes(y = AUC), position = position_dodge(0.9), vjust = -0.5)+
  guides(color='none')+coord_cartesian(ylim = c(-0.15, 1))+
  scale_color_manual(values=colormap)+
  theme(legend.position = c(0.75,0.25))+ 
  geom_line(size = 1)
  # geom_smooth(aes(x = FPR, y = TPR,colour = Method),size = 1,se = FALSE)

library(cowplot)
gg2 <- ggdraw() + 
  draw_plot(p4, 0,0,1, 1) + draw_plot(p3, 0.37,0.13,0.6,0.4)
print(gg2)


df = data.frame(weights = as.vector(Ksimilarity))
h4 = ggplot(df,aes(x = weights,y=..count../sum(..count..))) +
  geom_histogram(aes(y=..count../sum(..count..)),bins = 30,colour="black", fill="white")+
  geom_density(aes(y = ..density../30),alpha=.2, fill="#FF6666")+
  geom_vline(aes(xintercept=thresh4), colour="#990000", size = 1)+
  geom_vline(aes(xintercept=threshK), colour="#990000", linetype="dashed",size = 1)+
  labs(x='Weights',y="Density")+ggtitle("Krandom") +
  # geom_density(alpha=.2, fill="#FF6666") +
  theme(axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 12), axis.text.y = element_text(size = 12),
        title=element_text(size=12),legend.text = element_text(size = 12),
        plot.title = element_text(size = 12))

library(cowplot)
gg2 <- ggdraw() +
  draw_plot(gg1, 0,0,0.32, 1) + draw_plot(gg2, 0.33,0,0.32,1) + draw_plot(h4, 0.67,0,0.32,1)+
  draw_plot_label(c("A", "B","C"), c(0, 0.33,0.66), c(1,1,1), size = 15,colour = "black")
print(gg2)

#########################################################################


# 画图PR ROC曲线
library(ggplot2)
ggdatapr = data.frame(Recall = prcurve.K@x.values[[1]],
                      Precision = prcurve.K@y.values[[1]])
label = data.frame(
  Recall = 0.25,
  Precision = 0.75,
  label = paste('PR AUC =',round(prauc.K,2))
)
p1 = ggplot(data = ggdatapr, mapping = aes(x = Recall, y = Precision))+
  geom_line(size = 0.5,color = "black")+
  geom_ribbon(aes(ymin = 0.722, ymax = Precision), fill = "#f57670", alpha = 0.23)+
  geom_label(data = label, aes(label = label))+ggtitle("PR Curve") +
  theme(axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 12), axis.text.y = element_text(size = 12),
        title=element_text(size=12),legend.text = element_text(size = 12),
        plot.title = element_text(size = 12))+
  # geom_text(aes(y = AUC), position = position_dodge(0.9), vjust = -0.5)+
  # geom_smooth(size = 1)+#guides(color='none')+
  scale_color_manual(values=c("#DB7093","#6A5ACD","#6495ED",'#32CD32'))+
  theme(legend.position = c(0.55,0.15))

ggdataroc = data.frame(FPR = roccurve.K@x.values[[1]],
                       TPR = roccurve.K@y.values[[1]])
label = data.frame(
  FPR = 0.8,
  TPR = 0.1,
  label = paste('ROC AUC =',round(rocauc.K,2))
)
p2 = ggplot(data = ggdataroc, mapping = aes(x = FPR, y = TPR))+
  geom_line(size = 0.5,color = "black")+
  geom_ribbon(aes(ymin = 0, ymax = TPR), fill = "#f57670", alpha = 0.22)+
  geom_label(data = label, aes(label = label))+ggtitle("ROC Curve") +
  theme(axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 12), axis.text.y = element_text(size = 12),
        title=element_text(size=12),legend.text = element_text(size = 12),
        plot.title = element_text(size = 12))+
  # geom_text(aes(y = AUC), position = position_dodge(0.9), vjust = -0.5)+
  # geom_smooth(size = 1)+#guides(color='none')+
  scale_color_manual(values=c("#DB7093","#6A5ACD","#6495ED",'#32CD32'))+
  theme(legend.position = c(0.55,0.15))

df = data.frame(weights = as.vector(Ksimilarity))
h4 = ggplot(df,aes(x = weights,y=..count../sum(..count..))) +
  geom_histogram(aes(y=..count../sum(..count..)),bins = 30,colour="black", fill="white")+
  geom_density(aes(y = ..density../30),alpha=.2, fill="#FF6666")+
  geom_vline(aes(xintercept=thresh4), colour="#990000", size = 1)+
  geom_vline(aes(xintercept=threshK), colour="#990000", linetype="dashed",size = 1)+
  labs(y="density")+ggtitle("Krandom") +
  # geom_density(alpha=.2, fill="#FF6666") +
  theme(axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 12), axis.text.y = element_text(size = 12),
        title=element_text(size=12),legend.text = element_text(size = 12),
        plot.title = element_text(size = 12))

library(cowplot)
gg2 <- ggdraw() +
  draw_plot(p1, 0,0,0.32, 1) + draw_plot(p2, 0.33,0,0.32,1) + draw_plot(h4, 0.67,0,0.32,1)+
  draw_plot_label(c("A", "B","C"), c(0, 0.33,0.66), c(1,1,1), size = 15,colour = "black")
print(gg2)


######################################################plot net
library(GGally)
library(sna)
g = as.network.matrix(gnet,matrix.type ="adjacency",directed = FALSE)
g %v% 'name' = c(1:Nt)

xg = gplot.layout.kamadakawai(g, NULL)    #获取位置
g %v% "x" = xg[, 1]
g %v% "y" = xg[, 2]

p1 = ggnet2(g,vjust = -0.6,mode = c("x", "y"),color = "black",
       size = 1, edge.size = 1,edge.alpha = 1)+
  guides(color = 'none',size = 'none')+
  # geom_point( size = 12, color = "white") +
  geom_point( size = 12, alpha = 0.5) +
  geom_point( size = 9) +
  geom_text(aes(label = c('A21','A30','A31','A33','A40','A46','A49','A63','K03')), color = "white", fontface = "bold") +
  theme(title=element_text(size=15),legend.text = element_text(size = 15))

temp = KThreshmat
temp[which(temp!=0)] = temp[which(temp!=0)]-threshK
temp = temp*3
g2 = as.network.matrix(temp,matrix.type ="adjacency",directed = FALSE,ignore.eval = FALSE,names.eval = 'weights')
g2 %v% 'name' = c(1:Nt)

g2 %v% "x" = xg[, 1]
g2 %v% "y" = xg[, 2]

p3 = ggnet2(g2,vjust = -0.6,mode = c("x", "y"),color = "black",
       size = 1, edge.size = "weights",edge.alpha = 1)+
  guides(color = 'none',size = 'none')+
  # geom_point( size = 12, color = "white") +
  geom_point( size = 12, alpha = 0.5) +
  geom_point( size = 9) +
  geom_text(aes(label = c('A21','A30','A31','A33','A40','A46','A49','A63','K03')), color = "white", fontface = "bold") +
  theme(title=element_text(size=15),legend.text = element_text(size = 15))

library(cowplot)
hhh <- ggdraw() +
  draw_plot(p1, 0,0,0.33, 1) + draw_plot(p2, 0.34,0,0.33,1) + draw_plot(p3, 0.67,0,0.33,1) +
  draw_plot_label(c("A", "B","C"), c(0, 0.33,0.66), c(1,1,1), size = 15,colour = "black")
print(hhh)
