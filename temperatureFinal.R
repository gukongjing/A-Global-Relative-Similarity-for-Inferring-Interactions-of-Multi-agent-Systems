rm(list = ls())

library(hrbrthemes)

library(forcats)

require(sp) # the classes and methods that make up spatial ops in R
require(maptools) # tools for reading and manipulating spatial objects
require(rgeos)
require(rgdal)
require(gpclib)
require(scales)
require(dplyr)
gpclibPermit()

library(ncdf4)    #读nc数据
ncfile<-nc_open("data\\temperature\\X113.218.6.72.75.20.19.58.nc")   #air.day.1981-2010.ltm.nc
# names(ncfile$var)
t <- ncvar_get(ncfile, "time")
t0 = c(1:length(t))    #原始大小

longi <- ncvar_get(ncfile, "lon")
lati <- ncvar_get(ncfile, "lat", verbose = F)
temperature <- ncvar_get(ncfile, "air") # store the data in a 3-dimensional array

#NA value
# fillvalue <- ncatt_get(ncfile, "air", "_FillValue")
nc_close(ncfile) 
# temperature[temperature == fillvalue$value] <- NA

# 读地理数据
library(rgdal)
library(plyr)
library(sf)
American_map = readOGR("data/USA_map/STATES.SHP")
AD1 <- American_map@data
AD2 <- data.frame(id=rownames(AD1),AD1)
American_map1 <- fortify(American_map)    #提取美国各州轮廓
American_map_data <- plyr::join(American_map1,AD2, type = "full")   #合并轮廓和其他数据
American_map_data<-American_map_data[,1:12]     #去掉无关数据
American_data <-subset(American_map_data,STATE_NAME!='Alaska'& STATE_NAME!='Hawaii'
                       & STATE_NAME != "District of Columbia")   #去掉非本土州

#求各州的中点
centres1 = coordinates(American_map)
centres1 = centres1[-c(27,50,51),]
centres = data.frame(id = c(1:48), lon = centres1[,1]+360, lat = centres1[,2])

#各州中点的温度代表该州
piece = sapply(c(1:48),function(id){
  clongi = centres[id,2]
  index1 = ceiling((clongi-longi[1])/(longi[2]-longi[1]))
  clati = centres[id,3]
  index2 = ceiling((clati-lati[1])/(lati[2]-lati[1]))
  temperature[index1,index2,]
})

# ind = t(sapply(c(1:48),function(id){
#   clongi = centres[id,2]
#   index1 = floor((clongi-longi[1])/(longi[2]-longi[1]))
#   clati = centres[id,3]
#   index2 = floor((clati-lati[1])/(lati[2]-lati[1]))+1
#   c(index1,index2)
# }))

C = lapply(t0,function(day){
  # clustering on each piece
  maxl = length(unique(piece[day,]))
  # kk = sample(2:(maxl-1),1)
  kk=2
  clus = as.numeric(fct_inorder(as.factor(kmeans(piece[day,],centers = kk)$cluster)))
  D = matrix(0,nrow = 48, ncol = 48)             # initialize before each for
  for (ki in 1:kk) {
    index = which(clus == ki)
    D[index,index] = 1           # count for each piece of network
  }
  D
})

C0 = Reduce("+", C)
diag(C0) = 0
C = C0
# write.table (C0, file ="data\\temperature\\Csimple2.txt", sep =" ", row.names =FALSE, col.names =TRUE)

C = read.table('data\\temperature\\Csimple2.txt', sep =" ", header =TRUE)
C = as.matrix(C)

#####截断C
library(tidyverse) # 加载tidyverse包
library(mclust) # 加载mclust包

Cstand = C/max(C)
C0 = as.vector(Cstand)

histdata = data.frame(Weight = C0)

mod5 <- Mclust(C0,G=2)
plot(mod5, what = "classification")
summary(mod5)

C0 = rbind(C0,mod5[["classification"]])
C1 = C0[,order(C0[1,])]
C1 = C1[,which(C1[1,]<=mean(C1[1,which(C1[2,]==2)]))]   #去掉比第二个均值大的数
C1 = C1[,which(C1[1,]>=mean(C1[1,which(C1[2,]==1)]))]   #去掉比第一个均值小的数
thresh = (max(C1[1,which(C1[2,]==1)])+min(C1[1,which(C1[2,]==2)]))/2    #属于类别1的最大值与属于类别2的最小值取平均

# 设置..density..频率分布直方图
p1 = ggplot(histdata,aes(x = Weight, y =..count../length(C0))) +
  geom_histogram(aes(x = Weight, y =..count../length(C0)),stat="bin",binwidth=0.1,
                 fill= '#DDDDDD', color = 'black',size = 0.3)+
  theme(axis.title.x = element_text(size = 12,family = 'sans'), axis.text.x = element_text(angle = 45, hjust=1, 
                                                                                           size = 12,family = 'sans'),
        axis.title.y = element_text(size = 12,family = 'sans'), axis.text.y = element_text(size = 12,family = 'sans'),
        title=element_text(size=12,family = 'sans'))+
  labs(y='Density')+geom_vline(xintercept=thresh,color = '#990000',size =1)

d = Cstand
d[which(d<thresh)] = 0

library(corrplot)
corrplot(d)

################画图验证
C = d
library(RColorBrewer)
library(hrbrthemes)
library(ggtext)

library(forcats)
library(pacman)
p_load(assertthat,tidyverse,ggraph,igraph)


edges <- map_dfr(centres$id, function(id){
  to = which(C[id,] != 0)
  weight <- C[id,to]
  tibble(from=id, to=to, weight=weight)
})
# edges <- edges%>%mutate(category=as.factor(category))
g <- graph_from_data_frame(edges, directed = FALSE, vertices = centres)
edges_for_plot <- edges%>%
  inner_join(centres%>%dplyr::select(id, lon, lat),by=c("from"="id"))%>%
  dplyr::rename(x=lon, y=lat)%>%
  inner_join(centres%>%dplyr::select(id,lon,lat),by=c("to"="id"))%>%
  dplyr::rename(xend=lon,yend=lat)
assert_that(nrow(edges_for_plot)==nrow(edges))


my_colormap <- colorRampPalette(rev(brewer.pal(11,'Spectral')))(32)


library(raster)     #读nc数据
nc01 <- "data\\temperature\\X113.218.6.72.75.20.19.58.nc"
dset01 <- raster(nc01)
dset01_df <- as.data.frame(dset01,xy = TRUE)   ### as.vector(slice) 替换framework的最后一列
# head(dset01_df)
tem = sapply(t0, function(ttt){
  as.vector(temperature[, , ttt])
})
tem = apply(tem, 1, mean)
dset01_df$mean.Daily.Air.temperature.at.2.m = tem

centres$temperature = apply(piece,2,mean)

## https://zhuanlan.zhihu.com/p/88921255
# p2 = ggplot(centres)+
#   geom_tile(data = dset01_df, aes(x=x, y=y, fill=mean.Daily.Air.temperature.at.2.m),alpha=1)+
#   geom_polygon(data=American_data,aes(x=long+360,y=lat,group=group),colour="#293d6b",fill="white",alpha=.01,size=0.3)+
#   geom_curve(aes(x=x,y=y,xend=xend,yend=yend,color=weight),   #color=category,
#              data=edges_for_plot,size = 0.6,curvature = 0,alpha=1)+
#   scale_size_continuous(guide = 'none',range = c(0.25,2))+  # scale for edge widths
#   geom_point(aes(x=lon,y=lat,fill=temperature),  # draw nodes    ,size=weight
#              shape=21,color='black',stroke=0.5)+
#   scale_colour_gradient(low = "#696969", high = "black",limits=c(0,1),breaks=c(0,0.5,1))+
#   scale_fill_gradientn(colours = my_colormap,name="degK",breaks = c(275,285,295)) +
#   theme(panel.grid = element_blank(),
#         panel.background = element_blank(),
#         axis.text = element_blank(),
#         axis.ticks = element_blank(),
#         axis.title = element_blank(),
#         plot.title = element_markdown(hjust = 0.5,vjust = .5,color = "black",
#                                       size = 12, margin = margin(t = 1, b = 12)),
#         plot.subtitle = element_markdown(hjust = 0,vjust = .5,size=15),
#         plot.caption = element_markdown(face = 'bold',size = 12),
#         legend.position = 'bottom',
#         legend.title = element_text(size=12),
#         legend.text = element_text(size=12)
#   )+ggtitle("K=2")
p2 = ggplot(centres)+
  geom_tile(data = dset01_df, aes(x=x, y=y, fill=mean.Daily.Air.temperature.at.2.m),alpha=1)+
  geom_polygon(data=American_data,aes(x=long+360,y=lat,group=group),colour="#293d6b",fill="white",alpha=.01,size=0.3)+
  geom_curve(aes(x=x,y=y,xend=xend,yend=yend),   #color=category,
             data=edges_for_plot,size = 0.6,curvature = 0,alpha=1,color='black')+
  scale_size_continuous(guide = 'none',range = c(0.25,2))+  # scale for edge widths
  geom_point(aes(x=lon,y=lat,fill=temperature),  # draw nodes    ,size=weight
             shape=21,color='black',stroke=0.5)+
  scale_fill_gradientn(colours = my_colormap,name="degK",breaks = c(275,285,295)) +
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.title = element_markdown(hjust = 0.5,vjust = .5,color = "black",
                                      size = 12, margin = margin(t = 1, b = 12)),
        plot.subtitle = element_markdown(hjust = 0,vjust = .5,size=15),
        plot.caption = element_markdown(face = 'bold',size = 12),
        legend.position = 'bottom',
        legend.title = element_text(size=12),
        legend.text = element_text(size=12)
  )+ggtitle("K=2")
p2
###################################################################################################

library(hrbrthemes)

library(forcats)

require(sp) # the classes and methods that make up spatial ops in R
require(maptools) # tools for reading and manipulating spatial objects
require(rgeos)
require(rgdal)
require(gpclib)
require(scales)
require(dplyr)
gpclibPermit()

library(ncdf4)    #读nc数据
ncfile<-nc_open("data\\temperature\\X113.218.6.72.75.20.19.58.nc")   #air.day.1981-2010.ltm.nc
# names(ncfile$var)
t <- ncvar_get(ncfile, "time")
t0 = c(1:length(t))    #原始大小

longi <- ncvar_get(ncfile, "lon")
lati <- ncvar_get(ncfile, "lat", verbose = F)
temperature <- ncvar_get(ncfile, "air") # store the data in a 3-dimensional array

#NA value
# fillvalue <- ncatt_get(ncfile, "air", "_FillValue")
nc_close(ncfile) 
# temperature[temperature == fillvalue$value] <- NA

# 读地理数据
library(rgdal)
library(plyr)
library(sf)
American_map = readOGR("data/USA_map/STATES.SHP")
AD1 <- American_map@data
AD2 <- data.frame(id=rownames(AD1),AD1)
American_map1 <- fortify(American_map)    #提取美国各州轮廓
American_map_data <- plyr::join(American_map1,AD2, type = "full")   #合并轮廓和其他数据
American_map_data<-American_map_data[,1:12]     #去掉无关数据
American_data <-subset(American_map_data,STATE_NAME!='Alaska'& STATE_NAME!='Hawaii'
                       & STATE_NAME != "District of Columbia")   #去掉非本土州

#求各州的中点
centres1 = coordinates(American_map)
centres1 = centres1[-c(27,50,51),]
centres = data.frame(id = c(1:48), lon = centres1[,1]+360, lat = centres1[,2])

#各州中点的温度代表该州
piece = sapply(c(1:48),function(id){
  clongi = centres[id,2]
  index1 = ceiling((clongi-longi[1])/(longi[2]-longi[1]))
  clati = centres[id,3]
  index2 = ceiling((clati-lati[1])/(lati[2]-lati[1]))
  temperature[index1,index2,]
})

# ind = t(sapply(c(1:48),function(id){
#   clongi = centres[id,2]
#   index1 = floor((clongi-longi[1])/(longi[2]-longi[1]))
#   clati = centres[id,3]
#   index2 = floor((clati-lati[1])/(lati[2]-lati[1]))+1
#   c(index1,index2)
# }))

C = lapply(t0,function(day){
  # clustering on each piece
  maxl = length(unique(piece[day,]))
  # kk = sample(2:(maxl-1),1)
  kk=5
  clus = as.numeric(fct_inorder(as.factor(kmeans(piece[day,],centers = kk)$cluster)))
  D = matrix(0,nrow = 48, ncol = 48)             # initialize before each for
  for (ki in 1:kk) {
    index = which(clus == ki)
    D[index,index] = 1           # count for each piece of network
  }
  D
})

C0 = Reduce("+", C)
diag(C0) = 0
C = C0
# write.table (C0, file ="data\\temperature\\Csimple5.txt", sep =" ", row.names =FALSE, col.names =TRUE)

C = read.table('data\\temperature\\Csimple5.txt', sep =" ", header =TRUE)
C = as.matrix(C)

#####截断C
library(tidyverse) # 加载tidyverse包
library(mclust) # 加载mclust包

Cstand = C/max(C)
C0 = as.vector(Cstand)

histdata = data.frame(Weight = C0)
# ggplot(histdata,aes(x = x, y =..count..)) +
#   geom_histogram(aes(x = x, y =..count..),stat="bin",binwidth=0.1, boundary = 0)+
#   geom_text(aes(label=as.character(round(..count..,2))),stat="bin",binwidth=1,boundary = 0,vjust=-0.5)

mod5 <- Mclust(C0,G=2)
plot(mod5, what = "classification")
summary(mod5)

C0 = rbind(C0,mod5[["classification"]])
C1 = C0[,order(C0[1,])]
C1 = C1[,which(C1[1,]<=mean(C1[1,which(C1[2,]==2)]))]   #去掉比第二个均值大的数
C1 = C1[,which(C1[1,]>=mean(C1[1,which(C1[2,]==1)]))]   #去掉比第一个均值小的数
thresh = (max(C1[1,which(C1[2,]==1)])+min(C1[1,which(C1[2,]==2)]))/2    #属于类别1的最大值与属于类别2的最小值取平均

# 设置..density..频率分布直方图
p3 = ggplot(histdata,aes(x = Weight, y =..count../length(C0))) +
  geom_histogram(aes(x = Weight, y =..count../length(C0)),stat="bin",binwidth=0.1,
                 fill= '#DDDDDD', color = 'black',size = 0.3)+
  theme(axis.title.x = element_text(size = 12,family = 'sans'), axis.text.x = element_text(angle = 45, hjust=1, 
                                                                                           size = 12,family = 'sans'),
        axis.title.y = element_text(size = 12,family = 'sans'), axis.text.y = element_text(size = 12,family = 'sans'),
        title=element_text(size=12,family = 'sans'))+
  labs(y='Density')+geom_vline(xintercept=thresh,color = '#990000',size =1)

d = Cstand
d[which(d<thresh)] = 0

library(corrplot)
corrplot(d)

################画图验证
C = d
library(RColorBrewer)
library(hrbrthemes)
library(ggtext)

library(forcats)
library(pacman)
p_load(assertthat,tidyverse,ggraph,igraph)


edges <- map_dfr(centres$id, function(id){
  to = which(C[id,] != 0)
  weight <- C[id,to]
  tibble(from=id, to=to, weight=weight)
})
# edges <- edges%>%mutate(category=as.factor(category))
g <- graph_from_data_frame(edges, directed = FALSE, vertices = centres)
edges_for_plot <- edges%>%
  inner_join(centres%>%dplyr::select(id, lon, lat),by=c("from"="id"))%>%
  dplyr::rename(x=lon, y=lat)%>%
  inner_join(centres%>%dplyr::select(id,lon,lat),by=c("to"="id"))%>%
  dplyr::rename(xend=lon,yend=lat)
assert_that(nrow(edges_for_plot)==nrow(edges))


my_colormap <- colorRampPalette(rev(brewer.pal(11,'Spectral')))(32)


library(raster)     #读nc数据
nc01 <- "data\\temperature\\X113.218.6.72.75.20.19.58.nc"
dset01 <- raster(nc01)
dset01_df <- as.data.frame(dset01,xy = TRUE)   ### as.vector(slice) 替换framework的最后一列
# head(dset01_df)
tem = sapply(t0, function(ttt){
  as.vector(temperature[, , ttt])
})
tem = apply(tem, 1, mean)
dset01_df$mean.Daily.Air.temperature.at.2.m = tem

centres$temperature = apply(piece,2,mean)

## https://zhuanlan.zhihu.com/p/88921255
p4 = ggplot(centres)+
  geom_tile(data = dset01_df, aes(x=x, y=y, fill=mean.Daily.Air.temperature.at.2.m))+
  geom_polygon(data=American_data,aes(x=long+360,y=lat,group=group),colour="#293d6b",fill="white",alpha=.01,size=0.3)+
  geom_curve(aes(x=x,y=y,xend=xend,yend=yend),   #color=category,
             data=edges_for_plot,size = 0.6,curvature = 0,alpha=1,color='black')+
  scale_size_continuous(guide = 'none',range = c(0.25,2))+  # scale for edge widths
  geom_point(aes(x=lon,y=lat,fill=temperature),  # draw nodes    ,size=weight
             shape=21,color='black',stroke=0.5)+
  scale_fill_gradientn(colours = my_colormap,name="degK",breaks = c(275,285,295)) +
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.title = element_markdown(hjust = 0.5,vjust = .5,color = "black",
                                      size = 12, margin = margin(t = 1, b = 12)),
        plot.subtitle = element_markdown(hjust = 0,vjust = .5,size=15),
        plot.caption = element_markdown(face = 'bold',size = 12),
        legend.position = 'bottom',
        legend.title = element_text(size=12),
        legend.text = element_text(size=12)
  )+ggtitle("K=5")
p4
####################################################################################################################################

library(hrbrthemes)

library(forcats)

require(sp) # the classes and methods that make up spatial ops in R
require(maptools) # tools for reading and manipulating spatial objects
require(rgeos)
require(rgdal)
require(gpclib)
require(scales)
require(dplyr)
gpclibPermit()

library(ncdf4)    #读nc数据
ncfile<-nc_open("data\\temperature\\X113.218.6.72.75.20.19.58.nc")   #air.day.1981-2010.ltm.nc
# names(ncfile$var)
t <- ncvar_get(ncfile, "time")
t0 = c(1:length(t))    #原始大小

longi <- ncvar_get(ncfile, "lon")
lati <- ncvar_get(ncfile, "lat", verbose = F)
temperature <- ncvar_get(ncfile, "air") # store the data in a 3-dimensional array

#NA value
# fillvalue <- ncatt_get(ncfile, "air", "_FillValue")
nc_close(ncfile) 
# temperature[temperature == fillvalue$value] <- NA

# 读地理数据
library(rgdal)
library(plyr)
library(sf)
American_map = readOGR("data/USA_map/STATES.SHP")
AD1 <- American_map@data
AD2 <- data.frame(id=rownames(AD1),AD1)
American_map1 <- fortify(American_map)    #提取美国各州轮廓
American_map_data <- plyr::join(American_map1,AD2, type = "full")   #合并轮廓和其他数据
American_map_data<-American_map_data[,1:12]     #去掉无关数据
American_data <-subset(American_map_data,STATE_NAME!='Alaska'& STATE_NAME!='Hawaii'
                       & STATE_NAME != "District of Columbia")   #去掉非本土州

#求各州的中点
centres1 = coordinates(American_map)
centres1 = centres1[-c(27,50,51),]
centres = data.frame(id = c(1:48), lon = centres1[,1]+360, lat = centres1[,2])

#各州中点的温度代表该州
piece = sapply(c(1:48),function(id){
  clongi = centres[id,2]
  index1 = ceiling((clongi-longi[1])/(longi[2]-longi[1]))
  clati = centres[id,3]
  index2 = ceiling((clati-lati[1])/(lati[2]-lati[1]))
  temperature[index1,index2,]
})

# ind = t(sapply(c(1:48),function(id){
#   clongi = centres[id,2]
#   index1 = floor((clongi-longi[1])/(longi[2]-longi[1]))
#   clati = centres[id,3]
#   index2 = floor((clati-lati[1])/(lati[2]-lati[1]))+1
#   c(index1,index2)
# }))

C = lapply(t0,function(day){
  # clustering on each piece
  maxl = length(unique(piece[day,]))
  # kk = sample(2:(maxl-1),1)
  kk=10
  clus = as.numeric(fct_inorder(as.factor(kmeans(piece[day,],centers = kk)$cluster)))
  D = matrix(0,nrow = 48, ncol = 48)             # initialize before each for
  for (ki in 1:kk) {
    index = which(clus == ki)
    D[index,index] = 1           # count for each piece of network
  }
  D
})

C0 = Reduce("+", C)
diag(C0) = 0
C = C0
# write.table (C0, file ="data\\temperature\\Csimple10.txt", sep =" ", row.names =FALSE, col.names =TRUE)

C = read.table('data\\temperature\\Csimple10.txt', sep =" ", header =TRUE)
C = as.matrix(C)

#####截断C
library(tidyverse) # 加载tidyverse包
library(mclust) # 加载mclust包

Cstand = C/max(C)
C0 = as.vector(Cstand)

histdata = data.frame(Weight = C0)
# ggplot(histdata,aes(x = x, y =..count..)) +
#   geom_histogram(aes(x = x, y =..count..),stat="bin",binwidth=0.1, boundary = 0)+
#   geom_text(aes(label=as.character(round(..count..,2))),stat="bin",binwidth=1,boundary = 0,vjust=-0.5)


mod5 <- Mclust(C0,G=2)
plot(mod5, what = "classification")
summary(mod5)

C0 = rbind(C0,mod5[["classification"]])
C1 = C0[,order(C0[1,])]
C1 = C1[,which(C1[1,]<=mean(C1[1,which(C1[2,]==2)]))]   #去掉比第二个均值大的数
C1 = C1[,which(C1[1,]>=mean(C1[1,which(C1[2,]==1)]))]   #去掉比第一个均值小的数
thresh = (max(C1[1,which(C1[2,]==1)])+min(C1[1,which(C1[2,]==2)]))/2    #属于类别1的最大值与属于类别2的最小值取平均

# 设置..density..频率分布直方图
p5 = ggplot(histdata,aes(x = Weight, y =..count../length(C0))) +
  geom_histogram(aes(x = Weight, y =..count../length(C0)),stat="bin",binwidth=0.1,
                 fill= '#DDDDDD', color = 'black',size = 0.3)+
  theme(axis.title.x = element_text(size = 12,family = 'sans'), axis.text.x = element_text(angle = 45, hjust=1, 
                                                                                           size = 12,family = 'sans'),
        axis.title.y = element_text(size = 12,family = 'sans'), axis.text.y = element_text(size = 12,family = 'sans'),
        title=element_text(size=12,family = 'sans'))+
  labs(y='Density')+geom_vline(xintercept=thresh,color = '#990000',size =1)

d = Cstand
d[which(d<thresh)] = 0

library(corrplot)
corrplot(d)

################画图验证
C = d
library(RColorBrewer)
library(hrbrthemes)
library(ggtext)

library(forcats)
library(pacman)
p_load(assertthat,tidyverse,ggraph,igraph)


edges <- map_dfr(centres$id, function(id){
  to = which(C[id,] != 0)
  weight <- C[id,to]
  tibble(from=id, to=to, weight=weight)
})
# edges <- edges%>%mutate(category=as.factor(category))
g <- graph_from_data_frame(edges, directed = FALSE, vertices = centres)
edges_for_plot <- edges%>%
  inner_join(centres%>%dplyr::select(id, lon, lat),by=c("from"="id"))%>%
  dplyr::rename(x=lon, y=lat)%>%
  inner_join(centres%>%dplyr::select(id,lon,lat),by=c("to"="id"))%>%
  dplyr::rename(xend=lon,yend=lat)
assert_that(nrow(edges_for_plot)==nrow(edges))


my_colormap <- colorRampPalette(rev(brewer.pal(11,'Spectral')))(32)


library(raster)     #读nc数据
nc01 <- "data\\temperature\\X113.218.6.72.75.20.19.58.nc"
dset01 <- raster(nc01)
dset01_df <- as.data.frame(dset01,xy = TRUE)   ### as.vector(slice) 替换framework的最后一列
# head(dset01_df)
tem = sapply(t0, function(ttt){
  as.vector(temperature[, , ttt])
})
tem = apply(tem, 1, mean)
dset01_df$mean.Daily.Air.temperature.at.2.m = tem

centres$temperature = apply(piece,2,mean)

## https://zhuanlan.zhihu.com/p/88921255
p6 = ggplot(centres)+
  geom_tile(data = dset01_df, aes(x=x, y=y, fill=mean.Daily.Air.temperature.at.2.m))+
  geom_polygon(data=American_data,aes(x=long+360,y=lat,group=group),colour="#293d6b",fill="white",alpha=.01,size=0.3)+
  geom_curve(aes(x=x,y=y,xend=xend,yend=yend),   #color=category,
             data=edges_for_plot,size = 0.6,curvature = 0,alpha=1,color='black')+
  scale_size_continuous(guide = 'none',range = c(0.25,2))+  # scale for edge widths
  geom_point(aes(x=lon,y=lat,fill=temperature),  # draw nodes    ,size=weight
             shape=21,color='black',stroke=0.5)+
  scale_fill_gradientn(colours = my_colormap,name="degK",breaks = c(275,285,295)) +
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.title = element_markdown(hjust = 0.5,vjust = .5,color = "black",
                                      size = 12, margin = margin(t = 1, b = 12)),
        plot.subtitle = element_markdown(hjust = 0,vjust = .5,size=15),
        plot.caption = element_markdown(face = 'bold',size = 12),
        legend.position = 'bottom',
        legend.title = element_text(size=12),
        legend.text = element_text(size=12)
  )+ggtitle("K=10")

##################################################################################################################################

library(hrbrthemes)

library(forcats)

require(sp) # the classes and methods that make up spatial ops in R
require(maptools) # tools for reading and manipulating spatial objects
require(rgeos)
require(rgdal)
require(gpclib)
require(scales)
require(dplyr)
gpclibPermit()

library(ncdf4)    #读nc数据
ncfile<-nc_open("data\\temperature\\X113.218.6.72.75.20.19.58.nc")   #air.day.1981-2010.ltm.nc
# names(ncfile$var)
t <- ncvar_get(ncfile, "time")
t0 = c(1:length(t))    #原始大小

longi <- ncvar_get(ncfile, "lon")
lati <- ncvar_get(ncfile, "lat", verbose = F)
temperature <- ncvar_get(ncfile, "air") # store the data in a 3-dimensional array

#NA value
# fillvalue <- ncatt_get(ncfile, "air", "_FillValue")
nc_close(ncfile) 
# temperature[temperature == fillvalue$value] <- NA

# 读地理数据
library(rgdal)
library(plyr)
library(sf)
American_map = readOGR("data/USA_map/STATES.SHP")
AD1 <- American_map@data
AD2 <- data.frame(id=rownames(AD1),AD1)
American_map1 <- fortify(American_map)    #提取美国各州轮廓
American_map_data <- plyr::join(American_map1,AD2, type = "full")   #合并轮廓和其他数据
American_map_data<-American_map_data[,1:12]     #去掉无关数据
American_data <-subset(American_map_data,STATE_NAME!='Alaska'& STATE_NAME!='Hawaii'
                       & STATE_NAME != "District of Columbia")   #去掉非本土州

#求各州的中点
centres1 = coordinates(American_map)
centres1 = centres1[-c(27,50,51),]
centres = data.frame(id = c(1:48), lon = centres1[,1]+360, lat = centres1[,2])

#各州中点的温度代表该州
piece = sapply(c(1:48),function(id){
  clongi = centres[id,2]
  index1 = ceiling((clongi-longi[1])/(longi[2]-longi[1]))
  clati = centres[id,3]
  index2 = ceiling((clati-lati[1])/(lati[2]-lati[1]))
  temperature[index1,index2,]
})

# ind = t(sapply(c(1:48),function(id){
#   clongi = centres[id,2]
#   index1 = floor((clongi-longi[1])/(longi[2]-longi[1]))
#   clati = centres[id,3]
#   index2 = floor((clati-lati[1])/(lati[2]-lati[1]))+1
#   c(index1,index2)
# }))

C = lapply(t0,function(day){
  # clustering on each piece
  maxl = length(unique(piece[day,]))
  # kk = sample(2:(maxl-1),1)
  kk=20
  clus = as.numeric(fct_inorder(as.factor(kmeans(piece[day,],centers = kk)$cluster)))
  D = matrix(0,nrow = 48, ncol = 48)             # initialize before each for
  for (ki in 1:kk) {
    index = which(clus == ki)
    D[index,index] = 1           # count for each piece of network
  }
  D
})

C0 = Reduce("+", C)
diag(C0) = 0
C = C0
# write.table (C0, file ="data\\temperature\\Csimple20.txt", sep =" ", row.names =FALSE, col.names =TRUE)

C = read.table('data\\temperature\\Csimple20.txt', sep =" ", header =TRUE)
C = as.matrix(C)

#####截断C
library(tidyverse) # 加载tidyverse包
library(mclust) # 加载mclust包

Cstand = C/max(C)
C0 = as.vector(Cstand)

histdata = data.frame(Weight = C0)
# ggplot(histdata,aes(x = x, y =..count..)) +
#   geom_histogram(aes(x = x, y =..count..),stat="bin",binwidth=0.1, boundary = 0)+
#   geom_text(aes(label=as.character(round(..count..,2))),stat="bin",binwidth=1,boundary = 0,vjust=-0.5)

mod5 <- Mclust(C0,G=2)
plot(mod5, what = "classification")
summary(mod5)

C0 = rbind(C0,mod5[["classification"]])
C1 = C0[,order(C0[1,])]
C1 = C1[,which(C1[1,]<=mean(C1[1,which(C1[2,]==2)]))]   #去掉比第二个均值大的数
C1 = C1[,which(C1[1,]>=mean(C1[1,which(C1[2,]==1)]))]   #去掉比第一个均值小的数
thresh = (max(C1[1,which(C1[2,]==1)])+min(C1[1,which(C1[2,]==2)]))/2    #属于类别1的最大值与属于类别2的最小值取平均

# 设置..density..频率分布直方图
p7 = ggplot(histdata,aes(x = Weight, y =..count../length(C0))) +
  geom_histogram(aes(x = Weight, y =..count../length(C0)),stat="bin",binwidth=0.1,
                 fill= '#DDDDDD', color = 'black',size = 0.3)+
  theme(axis.title.x = element_text(size = 12,family = 'sans'), axis.text.x = element_text(angle = 45, hjust=1, 
                                                                                           size = 12,family = 'sans'),
        axis.title.y = element_text(size = 12,family = 'sans'), axis.text.y = element_text(size = 12,family = 'sans'),
        title=element_text(size=12,family = 'sans'))+
  labs(y='Density')+geom_vline(xintercept=thresh,color = '#990000',size =1)

d = Cstand
d[which(d<thresh)] = 0

library(corrplot)
corrplot(d)

################画图验证
C = d
library(RColorBrewer)
library(hrbrthemes)
library(ggtext)

library(forcats)
library(pacman)
p_load(assertthat,tidyverse,ggraph,igraph)


edges <- map_dfr(centres$id, function(id){
  to = which(C[id,] != 0)
  weight <- C[id,to]
  tibble(from=id, to=to, weight=weight)
})
# edges <- edges%>%mutate(category=as.factor(category))
g <- graph_from_data_frame(edges, directed = FALSE, vertices = centres)
edges_for_plot <- edges%>%
  inner_join(centres%>%dplyr::select(id, lon, lat),by=c("from"="id"))%>%
  dplyr::rename(x=lon, y=lat)%>%
  inner_join(centres%>%dplyr::select(id,lon,lat),by=c("to"="id"))%>%
  dplyr::rename(xend=lon,yend=lat)
assert_that(nrow(edges_for_plot)==nrow(edges))


my_colormap <- colorRampPalette(rev(brewer.pal(11,'Spectral')))(32)


library(raster)     #读nc数据
nc01 <- "data\\temperature\\X113.218.6.72.75.20.19.58.nc"
dset01 <- raster(nc01)
dset01_df <- as.data.frame(dset01,xy = TRUE)   ### as.vector(slice) 替换framework的最后一列
# head(dset01_df)
tem = sapply(t0, function(ttt){
  as.vector(temperature[, , ttt])
})
tem = apply(tem, 1, mean)
dset01_df$mean.Daily.Air.temperature.at.2.m = tem

centres$temperature = apply(piece,2,mean)

## https://zhuanlan.zhihu.com/p/88921255
p8 = ggplot(centres)+
  geom_tile(data = dset01_df, aes(x=x, y=y, fill=mean.Daily.Air.temperature.at.2.m))+
  geom_polygon(data=American_data,aes(x=long+360,y=lat,group=group),colour="#293d6b",fill="white",alpha=.01,size=0.3)+
  geom_curve(aes(x=x,y=y,xend=xend,yend=yend),   #color=category,
             data=edges_for_plot,size = 0.6,curvature = 0,alpha=1,color='black')+
  scale_size_continuous(guide = 'none',range = c(0.25,2))+  # scale for edge widths
  geom_point(aes(x=lon,y=lat,fill=temperature),  # draw nodes    ,size=weight
             shape=21,color='black',stroke=0.5)+
  scale_fill_gradientn(colours = my_colormap,name="degK",breaks = c(275,285,295)) +
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.title = element_markdown(hjust = 0.5,vjust = .5,color = "black",
                                      size = 12, margin = margin(t = 1, b = 12)),
        plot.subtitle = element_markdown(hjust = 0,vjust = .5,size=15),
        plot.caption = element_markdown(face = 'bold',size = 12),
        legend.position = 'bottom',
        legend.title = element_text(size=12),
        legend.text = element_text(size=12)
  )+labs(title="K=20")
p8

##################################################################################################

library(cowplot)
gg1 <- ggdraw() +
  draw_plot(p2, 0,0.4,0.25, 0.6) + draw_plot(p1, 0,0,0.25,0.4)+
  draw_plot(p4, 0.25,0.4,0.25, 0.6) + draw_plot(p3, 0.25,0,0.25,0.4)+
  draw_plot(p6, 0.5,0.4,0.25, 0.6) + draw_plot(p5, 0.5,0,0.25,0.4)+
  draw_plot(p8, 0.75,0.4,0.25, 0.6) + draw_plot(p7, 0.75,0,0.25,0.4)
  
print(gg1)

