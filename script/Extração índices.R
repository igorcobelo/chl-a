#Packages

library(raster)
library(sf)
library(FIELDimageR)
library(ggplot2)


#Read raster
setwd("C:/Users/Higor Ferreira/Documents/VANT- ULTIMO/VANT- ULTIMO/20m/UEG_201119_20m/3_dsm_ortho/2_mosaic")
EX1<-raster::stack("ueg.tif")

#Verify spatial resolution
res(EX1)


#Plot raster
plotRGB(EX1, r = 1, g = 2, b= 3)



#Read shapefile of mesocosm
cxs <- sf::st_read("polygons_cx.shp")

#Plot shp
plot(cxs$geometry)

#Clip mesocosm in raster
rr <- mask(EX1, cxs)

#Extract values RGB
ex <- raster::extract(rr, cxs, sp=T)

#Create a data frame with the RGB mean values
df = as.data.frame(do.call(rbind, lapply(ex, colMeans))) 


####Extract vegetation indices

#NGRDI
NGRDI <- (rr[[2]]-rr[[1]])/(rr[[2]]+rr[[1]])
ex_NGRDI <- raster::extract(NGRDI, cxs, sp=T,fun=mean,df=T)
ex_NGRDI$layer

#NGBDI
NGBDI <- (rr[[2]]-rr[[3]])/(rr[[2]]+rr[[3]])
ex_NGBDI <- raster::extract(NGBDI, cxs, sp=T,fun=mean,df=T)
ex_NGBDI$layer

#GLI
GLI <- (2*rr[[2]]-rr[[1]]-rr[[3]])/(2*rr[[2]]+rr[[1]]+rr[[3]])
ex_GLI <- raster::extract(GLI, cxs, sp=T,fun=mean,df=T)
ex_GLI$layer

#EXG
EXG <- 2*rr[[2]]-rr[[1]]-rr[[3]]
ex_EXG <- raster::extract(EXG, cxs, sp=T,fun=mean,df=T)
ex_EXG$layer

#VARI
VARI <- (rr[[2]]-rr[[1]])/(rr[[2]]+rr[[1]]-rr[[3]])
ex_VARI <- raster::extract(VARI, cxs, sp=T,fun=mean,df=T)
ex_VARI$layer

#VWRI
VWRI <- (rr[[2]]-rr[[3]]-rr[[1]])/(rr[[1]]+rr[[3]]+rr[[2]])
ex_VWRI <- raster::extract(VWRI, cxs, sp=T,fun=mean,df=T)
ex_VWRI$layer

#SCI
SCI <- (rr[[1]]-rr[[2]])/(rr[[1]]+rr[[2]])
ex_SCI <- raster::extract(SCI, cxs, sp=T,fun=mean,df=T)
ex_SCI$layer

#SI
SI <- (rr[[1]]-rr[[3]])/(rr[[1]]+rr[[3]])
ex_SI <- raster::extract(SI, cxs, sp=T,fun=mean,df=T)
ex_SI$layer

#Create a data frame with indices values

df <- data.frame(Red=df$ueg.1, Green=df$ueg.2, Blue=df$ueg.3,
           EXG=ex_EXG$layer, GLI=ex_GLI$layer,NGBDI=ex_NGBDI$layer,
           NGRDI=ex_NGRDI$layer,SCI=ex_SCI$layer,SI=ex_SI$layer,
           VARI=ex_VARI$layer,VWRI=ex_VWRI$layer)
write.table(df, "novos_indices.csv",sep=';')


#####Crop and rotate images for each index (for plot)
#1) Crop two times (because the image is big)
#2) rotate with a specific rotation angle

NGRDI_crop <- fieldCrop(mosaic = NGRDI, fast.plot = T)
NGRDI_crop <- fieldCrop(mosaic = NGRDI_crop, fast.plot = T)
NGRDI_rotated<-fieldRotate(mosaic = NGRDI_crop, theta = -3.27,clockwise = T) # h=horizontal

NGBDI_crop <- fieldCrop(mosaic = NGBDI, fast.plot = T)
NGBDI_crop <- fieldCrop(mosaic = NGBDI_crop, fast.plot = T)
NGBDI_rotated<-fieldRotate(mosaic = NGBDI_crop, theta = -3.27,clockwise = T) # h=horizontal

GLI_crop <- fieldCrop(mosaic = GLI, fast.plot = T)
GLI_crop <- fieldCrop(mosaic = GLI_crop, fast.plot = T)
GLI_rotated<-fieldRotate(mosaic = GLI_crop, theta = -3.27,clockwise = T) # h=horizontal

EXG_crop <- fieldCrop(mosaic = EXG, fast.plot = T)
EXG_crop <- fieldCrop(mosaic = EXG_crop, fast.plot = T)
EXG_rotated<-fieldRotate(mosaic = EXG_crop, theta = -3.27,clockwise = T) # h=horizontal

VARI_crop <- fieldCrop(mosaic = VARI, fast.plot = T)
VARI_crop <- fieldCrop(mosaic = VARI_crop, fast.plot = T)
VARI_rotated<-fieldRotate(mosaic = VARI_crop, theta = -3.27,clockwise = T) # h=horizontal

VWRI_crop <- fieldCrop(mosaic = VWRI, fast.plot = T)
VWRI_crop <- fieldCrop(mosaic = VWRI_crop, fast.plot = T)
VWRI_rotated<-fieldRotate(mosaic = VWRI_crop, theta = -3.27,clockwise = T) # h=horizontal

SCI_crop <- fieldCrop(mosaic = SCI, fast.plot = T)
SCI_crop <- fieldCrop(mosaic = SCI_crop, fast.plot = T)
SCI_rotated<-fieldRotate(mosaic = SCI_crop, theta = -3.27,clockwise = T) # h=horizontal


SI_crop <- fieldCrop(mosaic = SI, fast.plot = T)
SI_crop <- fieldCrop(mosaic = SI_crop, fast.plot = T)
SI_rotated<-fieldRotate(mosaic = SI_crop, theta = -3.27,clockwise = T) # h=horizontal


########Plots



#NGBDI

# convert to a df for plotting in two steps,
# First, to a SpatialPointsDataFrame
NGBDI_rotated_pts <- rasterToPoints(NGBDI_rotated, spatial = TRUE)
# Then to a 'conventional' dataframe
NGBDI_rotated_df  <- data.frame(NGBDI_rotated_pts)
colnames(NGBDI_rotated_df)[1] <- "NGBDI"

ngbdi_plot <- ggplot() +
  geom_tile(data = NGBDI_rotated_df , aes(x = x, y = y, fill = NGBDI)) + 
  scale_fill_gradientn(colours=c("red","yellow","forestgreen"))+
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_blank(),axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),axis.title.x = element_blank(),
        axis.title.y = element_blank())
  

#NGRDI

# convert to a df for plotting in two steps,
# First, to a SpatialPointsDataFrame
NGRDI_rotated_pts <- rasterToPoints(NGRDI_rotated, spatial = TRUE)
# Then to a 'conventional' dataframe
NGRDI_rotated_df  <- data.frame(NGRDI_rotated_pts)
colnames(NGRDI_rotated_df)[1] <- "NGRDI"

ngrdi_plot <- ggplot() +
  geom_tile(data = NGRDI_rotated_df , aes(x = x, y = y, fill = NGRDI)) + 
  scale_fill_gradientn(colours=c("red","yellow","forestgreen"))+
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_blank(),axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),axis.title.x = element_blank(),
        axis.title.y = element_blank())

#GLI

# convert to a df for plotting in two steps,
# First, to a SpatialPointsDataFrame
GLI_rotated_pts <- rasterToPoints(GLI_rotated, spatial = TRUE)
# Then to a 'conventional' dataframe
GLI_rotated_df  <- data.frame(GLI_rotated_pts)
colnames(GLI_rotated_df)[1] <- "GLI"

gli_plot <- ggplot() +
  geom_tile(data = GLI_rotated_df , aes(x = x, y = y, fill = GLI)) + 
  scale_fill_gradientn(colours=c("red","yellow","forestgreen"))+
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_blank(),axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),axis.title.x = element_blank(),
        axis.title.y = element_blank())

#EXG

# convert to a df for plotting in two steps,
# First, to a SpatialPointsDataFrame
EXG_rotated_pts <- rasterToPoints(EXG_rotated, spatial = TRUE)
# Then to a 'conventional' dataframe
EXG_rotated_df  <- data.frame(EXG_rotated_pts)
colnames(EXG_rotated_df)[1] <- "EXG"

EXG_plot <- ggplot() +
  geom_tile(data = EXG_rotated_df , aes(x = x, y = y, fill = EXG)) + 
  scale_fill_gradientn(colours=c("red","yellow","forestgreen"))+
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_blank(),axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),axis.title.x = element_blank(),
        axis.title.y = element_blank())

#VARI

# convert to a df for plotting in two steps,
# First, to a SpatialPointsDataFrame
VARI_rotated_pts <- rasterToPoints(VARI_rotated, spatial = TRUE)
# Then to a 'conventional' dataframe
VARI_rotated_df  <- data.frame(VARI_rotated_pts)
colnames(VARI_rotated_df)[1] <- "VARI"

VARI_plot <- ggplot() +
  geom_tile(data = VARI_rotated_df , aes(x = x, y = y, fill = VARI)) + 
  scale_fill_gradientn(colours=c("red","yellow","forestgreen"))+
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_blank(),axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),axis.title.x = element_blank(),
        axis.title.y = element_blank())


#VWRI

# convert to a df for plotting in two steps,
# First, to a SpatialPointsDataFrame
VWRI_rotated_pts <- rasterToPoints(VWRI_rotated, spatial = TRUE)
# Then to a 'conventional' dataframe
VWRI_rotated_df  <- data.frame(VWRI_rotated_pts)
colnames(VWRI_rotated_df)[1] <- "VWRI"

VWRI_plot <- ggplot() +
  geom_tile(data = VWRI_rotated_df , aes(x = x, y = y, fill = VWRI)) + 
  scale_fill_gradientn(colours=c("red","yellow","forestgreen"))+
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_blank(),axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),axis.title.x = element_blank(),
        axis.title.y = element_blank())


#SCI

# convert to a df for plotting in two steps,
# First, to a SpatialPointsDataFrame
SCI_rotated_pts <- rasterToPoints(SCI_rotated, spatial = TRUE)
# Then to a 'conventional' dataframe
SCI_rotated_df  <- data.frame(SCI_rotated_pts)
colnames(SCI_rotated_df)[1] <- "SCI"

SCI_plot <- ggplot() +
  geom_tile(data = SCI_rotated_df , aes(x = x, y = y, fill = SCI)) + 
  scale_fill_gradientn(colours=c("red","yellow","forestgreen"))+
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_blank(),axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),axis.title.x = element_blank(),
        axis.title.y = element_blank())


#SI

# convert to a df for plotting in two steps,
# First, to a SpatialPointsDataFrame
SI_rotated_pts <- rasterToPoints(SI_rotated, spatial = TRUE)
# Then to a 'conventional' dataframe
SI_rotated_df  <- data.frame(SI_rotated_pts)
colnames(SI_rotated_df)[1] <- "SI"

SI_plot <- ggplot() +
  geom_tile(data = SI_rotated_df , aes(x = x, y = y, fill = SI)) + 
  scale_fill_gradientn(colours=c("red","yellow","forestgreen"))+
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_blank(),axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),axis.title.x = element_blank(),
        axis.title.y = element_blank())