
#Load the relevant libraries.  
library(ggplot2)
library(purrr)
library(scrattch.vis)
library(scrattch.hicat)
library(Matrix)
library(dplyr)
library(tidyr)
library(parallel)
library(bigstatsr)
library(scrattch.bigcat)#, lib.loc = "/allen/programs/celltypes/workgroups/mct-t200/Cindy_analysis/my_Tools/myPY/envs/r4_env/lib/R/library")
library(data.table)
library(arrow)
library(BiocNeighbors)


indir = "/allen/programs/celltypes/workgroups/mct-t200/Cindy_analysis/my_Tools/TEST/script_examples/2023_constellation"
outdir = d ="/allen/programs/celltypes/workgroups/mct-t200/Cindy_analysis/my_Tools/cy_plotfolio/Centroid_app"

###############
# load data 
###############

load(file.path(indir,"LSX.anno.rda"))
load(file.path(indir,"LSX.cl.df.rda"))
source(file.path(indir,"plot_constellation.R"))

############################
# cluster constellation 
###############

# set cluster to cell vector
cl.cl = with(anno %>% filter(platform=="10X_cells_v3" & 
                               !is.na(subclass_id_label) & 
                               class_id_label == "07 LSX GABA"), 
             setNames(cluster_id, sample_id))

# load reduced dimension data (PCA in this case)
rd.dat.df = read_parquet("LSX.rd.parquet")
rd.dat = as.matrix(rd.dat.df[,-1])
row.names(rd.dat) = rd.dat.df[[1]]
rd.dat = rd.dat[intersect(rownames(rd.dat),names(cl.cl)),]

# calculate and summarize distance between all cells
result = get_knn_graph(rd.dat, cl=cl.cl, k =50)
save(result, file="cluster.knn.graph.result.rda")

load(file.path(indir,"cluster.knn.graph.result.rda"))
knn.cl.df = result$knn.cl.df
knn.cl.df = knn.cl.df %>% group_by(cl.from) %>% mutate(cl.from.rank = rank(-Freq))
knn.cl.df = knn.cl.df %>% group_by(cl.to) %>% mutate(cl.to.rank = rank(-Freq))
# criteria used to filter non-meaningful connections
## this depends a bit on the size and complexity of the dataset. Try a few parameters and check the plot to see which edges are robust and make sense and which are not.
select.knn.cl.df = with(knn.cl.df, knn.cl.df[ (frac > 0.03 | frac > 0.003 & Freq > 25) &  (cl.from.rank < 4| cl.to.rank < 4),])

select.knn.cl.df = with(knn.cl.df, knn.cl.df[odds > 1 & pval.log < log(1/100000) & (frac > 0.1 | frac > 0.03 & Freq > 100) & (cl.from.rank < 4| cl.to.rank < 4),])


prefix = "LSX"
umap.fn = file.path(indir,paste0(prefix, ".umap.2d.csv"))
umap.df = as.data.frame(fread(umap.fn,header=TRUE))
row.names(umap.df) = umap.df[[1]] 
umap.df = umap.df[,-1]


# make a dataframe where each row represents a node in the constellation. In this case cluster level. 
# Note1: to have the node size scaled to the size of the group the function is expecting "cluster_size"
# Note2: colors are included in the dataframe and picked up by the function 
cl.df <- anno %>% 
  group_by(cluster_id, cluster_label, cluster_color,
           supertype_id, supertype_id_label, supertype_color, 
           subclass_id, subclass_id_label, subclass_color,
           class_id,class_id_label,class_color) %>% 
  summarise(gene.counts = mean(gene.counts.0),
            umi.counts = mean(umi.counts),
            qc.score = mean(qc.score),
            cluster_size=n()
            ) 

# get centroid position of each group in umap space
cl.center.df = as.data.frame(get_RD_cl_center(umap.df, cl.cl))
cl.center.df$cl = row.names(cl.center.df)
cl.center.df$cluster_id <- cl.center.df$cl
cl.df$cluster_id <- as.character(cl.df$cluster_id)
# join all other group metadata
cl.center.df <- cl.center.df %>% left_join(cl.df)
#cl.center.df$cluster_color <- cl.center.df$supertype_color 

# select only clusters that you want to plot (in this case all clusters)
tmp.cl = cl.center.df$cl
tmp.knn.cl.df = select.knn.cl.df  %>% filter(cl.from %in% tmp.cl & cl.to %in% tmp.cl)

# ready to plot a constellation!
c.plot=plot_constellation(tmp.knn.cl.df, 
                          cl.center.df=cl.center.df, 
                          out.dir=outdir,
                          prefix=paste0(prefix,"_cluster_dodge"),
                          node.label="cluster_id",
                          exxageration=4,
                          label.size=4,
                          plot.parts=FALSE,
                          return.list = T,
                          label_repel = T,
                          node.dodge = T 
)

# lets change the color of the nodes to subclasses
cl.center.df$cluster_color <- cl.center.df$subclass_color

c.plot=plot_constellation(tmp.knn.cl.df, 
                          cl.center.df=cl.center.df, 
                          out.dir=outdir,
                          prefix=paste0(prefix,"_cluster_subclasscol"),
                          node.label="cluster_id",
                          exxageration=4,
                          label.size=4,
                          plot.parts=FALSE,
                          return.list = T,
                          label_repel = F,
                          node.dodge = T,
)

# there is an option to highlight nodes and/or highlight (or dim) specific edges
# to highlight nodes set parameters "highlight..."
# to manipulate edges set paremetes "edge_...". You can further adjust fg.alpha and/or bg.alpha of edges.
c.plot=plot_constellation(tmp.knn.cl.df, 
                          cl.center.df=cl.center.df, 
                          out.dir=outdir,
                          prefix=paste0(prefix,"_cluster_subclasscol_highlight"),
                          node.label="cluster_id",
                          exxageration=4,
                          label.size=4,
                          plot.parts=FALSE,
                          return.list = T,
                          label_repel = F,
                          node.dodge = T,
                          highlight_nodes = cl.center.df$cl[1:15],
                          highlight_color = "firebrick1",
                          highlight_width=2,
                          highlight_labsize = 10,
                          edge_mark_list = cl.center.df$cl[75:85],
                          edge_marking = "highlight"
                          
)

# another option for plotting is to plot a "hull" around a specific group. In this case supertype_id is used to plot a boundary around the group. This feature is the hardest to perfect but can be useful as a guide for exploring complex plots.

cl.center.df$clade_id <- cl.center.df$subclass_id
cl.center.df$clade_color <- cl.center.df$subclass_color

c.plot=plot_constellation(tmp.knn.cl.df, 
                          cl.center.df=cl.center.df, 
                          out.dir=outdir,
                          prefix=paste0(prefix,"_cluster_subclasscol_hull"),
                          node.label="cluster_id",
                          exxageration=4,
                          label.size=4,
                          plot.parts=FALSE,
                          return.list = T,
                          label_repel = F,
                          node.dodge = T,
                          plot.hull = cl.center.df$clade_id
)


############################
# supertype constellation 
###############

# ok one final example. Instead of plotting at a cluster level you can plot any other level so lets use supertype.

### supertype constellation for global dataset
st.cl = with(anno %>% filter(platform=="10X_cells_v3" & 
                               !is.na(subclass_id_label) & 
                               class_id_label == "07 LSX GABA"), 
             setNames(supertype_id, sample_id))

rd.dat.df = read_parquet("LSX.rd.parquet")
rd.dat = as.matrix(rd.dat.df[,-1])
row.names(rd.dat) = rd.dat.df[[1]]
rd.dat = rd.dat[intersect(rownames(rd.dat),names(st.cl)),]

result = get_knn_graph(rd.dat, cl=st.cl, k =50)
save(result, file="supertype.knn.graph.result.rda")

#load("supertype.knn.graph.result.rda")
knn.cl.df = result$knn.cl.df
knn.cl.df = knn.cl.df %>% group_by(cl.from) %>% mutate(cl.from.rank = rank(-Freq))
knn.cl.df = knn.cl.df %>% group_by(cl.to) %>% mutate(cl.to.rank = rank(-Freq))
select.knn.cl.df = with(knn.cl.df, knn.cl.df[ (frac > 0.03 | frac > 0.003 & Freq > 25) &  (cl.from.rank < 4| cl.to.rank < 4),])

prefix = "LSX"
umap.fn = file.path(d,paste0(prefix, ".umap.2d.csv"))
umap.df = as.data.frame(fread(umap.fn,header=TRUE))
row.names(umap.df) = umap.df[[1]] 
umap.df = umap.df[,-1]

super.df <- anno %>% 
  group_by( supertype_id, supertype_id_label, supertype_color, 
            subclass_id, subclass_id_label, subclass_color,
            class_id,class_id_label,class_color) %>% 
  summarise(cluster_size=n()) 

cl.center.df = as.data.frame(get_RD_cl_center(umap.df,st.cl))
cl.center.df$cl = row.names(cl.center.df)
cl.center.df$supertype_id <- cl.center.df$cl
super.df$supertype_id <- as.character(super.df$supertype_id)
cl.center.df <- cl.center.df %>% left_join(super.df)
cl.center.df$cluster_color <- cl.center.df$supertype_color 
cl.center.df$cluster_id <- cl.center.df$supertype_id
cl.center.df$cluster_label <- cl.center.df$supertype_id_label


tmp.cl = cl.center.df$cl
tmp.knn.cl.df = select.knn.cl.df  %>% filter(cl.from %in% tmp.cl & cl.to %in% tmp.cl)

c.plot=plot_constellation(tmp.knn.cl.df, 
                          cl.center.df=cl.center.df, 
                          out.dir=outdir,
                          prefix=paste0(prefix,"_supertype_dodge_alt"),
                          node.label="cluster_id",
                          exxageration=4,
                          label.size=4,
                          plot.parts=FALSE,
                          return.list = T,
                          label_repel = F,
                          node.dodge = T
)

cl.center.df$cluster_color <- cl.center.df$class_color

c.plot=plot_constellation(tmp.knn.cl.df, 
                          cl.center.df=cl.center.df, 
                          out.dir=outpath,
                          prefix=paste0(prefix,"_supertype_classcol"),
                          node.label="cluster_id",
                          exxageration=4,
                          label.size=4,
                          plot.parts=FALSE,
                          return.list = T,
                          label_repel = F,
                          node.dodge = T
)
