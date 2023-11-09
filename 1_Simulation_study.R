library(spatstat)
library(spatstat.linnet)
library(spatstat.geom)
library(rgeos)
library(spdep)
library(maptools)

# Set the working directory...
# setwd("")

source("Functions/sim_pattern_lpp.R")

load("Data/grafo_intersecciones.rda")
lines=as.SpatialLines.psp(as.psp(grafo_intersecciones))
segment_lengths=SpatialLinesLengths(lines)
proj4string(lines)="+proj=utm +zone=30 +ellps=WGS84 +units=m +no_defs"
lines=SpatialLinesDataFrame(lines,data.frame(Line=1:length(lines)),match.ID = F)
proj4string(lines)="+proj=utm +zone=30 +ellps=WGS84 +units=m +no_defs"

# Simulate data

set.seed(12345)

load("Data/X1.rda")
load("Data/X2.rda")
alpha=-4
beta1=0.5
beta2=1

for (error in c(0,5,10,20)){

  for (sim in 1:100){
    
    print(sim)
    
    # Simulate 'true' pattern
    
    if (!file.exists(paste0("SimData/sim_pattern_",sim,".rda"))){
      sim_pattern = sim_pattern_lpp(X1,X2,alpha,beta1,beta2,segment_lengths,linnet_object=grafo_intersecciones)
      save(sim_pattern,file=paste0("SimData/sim_pattern_",sim,".rda"))
    } else {
      load(paste0("SimData/sim_pattern_",sim,".rda"))
    }
    
    # Extract coordinates and add noise
    
    points=data.frame(x=sim_pattern[,1],y=sim_pattern[,2])
    coordinates(points)=~x+y
    proj4string(points)="+proj=utm +zone=30 +ellps=WGS84 +units=m +no_defs"
    # plot(points,pch=19,cex=0.5)
  
    points@coords[,1]=points@coords[,1]+rnorm(length(points),0,error)
    points@coords[,2]=points@coords[,2]+rnorm(length(points),0,error)
    # plot(points,pch=19,cex=0.5,col="red",add=T)
    save(points,file=paste0("SimData/points_moved_",error,"_",sim,".rda"))
    
    # Check network dimension
    
    grafo_intersecciones$lines$n
    
    # Compute network's centroids
    
    network_centroids=midpoints.psp(as.psp(grafo_intersecciones))
    
    # Compute distances between events and segments
    
    dist_events_segments=gDistance(lines,points,byid = T)
    
    # For each event, compute the 4 minimum distances and the segments corresponding to these distances
    
    distances=c()
    segments=c()
    weights=c()
    weights_squared=c()
    K=4 # choice
    for (i in 1:length(points)){
      distances_sort=sort(dist_events_segments[i,])
      segments_sort=as.numeric(names(distances_sort)[1:K])
      distances=rbind(distances,distances_sort[1:K])
      segments=rbind(segments,segments_sort)
      weights=rbind(weights,1/distances_sort[1:K]/(sum(1/distances_sort[1:K])))
      weights_squared=rbind(weights_squared,1/distances_sort[1:K]^2/(sum(1/distances_sort[1:K]^2)))
    }
    colnames(distances)=c("S1","S2","S3","S4")
    colnames(segments)=c("S1","S2","S3","S4")
    colnames(weights)=c("S1","S2","S3","S4")
    colnames(weights_squared)=c("S1","S2","S3","S4")
    
    save(distances,file=paste0("SimData/distances_",K,"_",error,"_sim",sim,".rda"))
    save(segments,file=paste0("SimData/segments_",K,"_",error,"_sim",sim,".rda"))
    save(weights,file=paste0("SimData/weights_",K,"_",error,"_sim",sim,".rda"))
    save(weights_squared,file=paste0("SimData/weights_squared_",K,"_",error,"_sim",sim,".rda"))
    
  }

}
