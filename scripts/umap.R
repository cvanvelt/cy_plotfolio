#library(scrattch.bigcat)
#library(data.table)

#load("anno.df.rda")
#load("rd.dat.list.rda")
#load("clean.sampled.cells.rda")
#load("cl.final.rda")
#load("rm.eigen.rda")
#rd.dat = filter_RD(rd.dat.list[["10X_cells"]]$rd.dat, rm.eigen, 0.7)

#dest.d="Manuscript/common/umap_constellation"
#for(n in c("CGE","MGE","PT","NP_L6CT_L6b","DG_SUB_CA")){


plot_umap_constellation <- function(umap.2d, cl, cl.df, select.knn.cl.df, dest.d=".", prefix="",...)
  {
    
    cl.center.df = as.data.frame(get_RD_cl_center(umap.2d,cl))
    cl.center.df$cl = row.names(cl.center.df)
    cl.center.df$cluster_id <- cl.df$cluster_id[match(cl.center.df$cl, cl.df$cl)]
    cl.center.df$cluster_color <- cl.df$cluster_color[match(cl.center.df$cl, cl.df$cl)]
    cl.center.df$cluster_label <- cl.df$cluster_label[match(cl.center.df$cl, cl.df$cl)] 
    cl.center.df$cluster_size <- cl.df$cluster_size[match(cl.center.df$cl, cl.df$cl)]
    tmp.cl = row.names(cl.center.df)
    tmp.knn.cl.df = select.knn.cl.df  %>% filter(cl.from %in% tmp.cl & cl.to %in% tmp.cl)
    p=plot_constellation(tmp.knn.cl.df, cl.center.df, node.label="cluster_id", out.dir=file.path(dest.d), prefix=prefix,...)    
    
    return(list(tmp.knn.cl.df, cl.center.df, p))
  }

find_umap <- function(select.cl, prefix, cl, anno.df=NULL, rd.dat, sampled.cells, dest.d, do.init=TRUE,meta.fields=NULL, nc=2, nn=25, md=0.3, niter=1)
  {
    tmp.sampled.cells= intersect(sampled.cells, names(cl[cl %in% select.cl]))
    rd.fn= paste0(dest.d,"/",prefix,".rd.dat.ref.10X.sampled.csv")
    fwrite(as.data.frame(rd.dat[tmp.sampled.cells,]), file=rd.fn, row.names=TRUE)
    
    umap.fn = paste0(dest.d,"/",prefix,".umap.",nc,"d.sampled.csv")  
    if(do.init){
      rd.cl.fn= paste0(dest.d,"/",prefix,".rd.dat.cl.means.csv")
      rd.cl.means = t(get_cl_means(rd.dat, cl))
      write.csv(rd.cl.means[select.cl,], file=rd.cl.fn)
  
      umap.cl.fn = paste0(dest.d,"/",prefix,".umap.",nc, "d.cl.means.csv")  
      #cmd = paste("~/zizhen/bin/run_umap.py -i", rd.cl.fn, "-o", umap.cl.fn, "-n 10")
      cmd = paste("/allen/programs/celltypes/workgroups/mct-t200/Cindy_analysis/my_Tools/run_umap.py -i", rd.cl.fn, "-o", umap.cl.fn, "-c", nc)
      print(cmd)
      system(cmd)
      
      umap.cl.fn = paste0(dest.d,"/",prefix,".umap.", nc, "d.cl.means.csv")  
      cl.center.df <- read.csv(umap.cl.fn)
      colnames(cl.center.df)=c("cl", paste0("Dim", 1:nc))
      cl.center.df$cl = as.character(cl.center.df$cl)
      umap.init = cl.center.df %>% left_join(data.frame(cl=cl, sample_id=names(cl)))
      row.names(umap.init)= umap.init$sample_id
      umap.init = umap.init[,paste0("Dim",1:nc)]
      umap.init.fn = paste0(dest.d,"/",prefix,".umap.init.csv")
      fwrite(umap.init[tmp.sampled.cells,], file= umap.init.fn, row.names=TRUE)
      cmd = paste("/allen/programs/celltypes/workgroups/mct-t200/Cindy_analysis/my_Tools/run_umap.py -i", rd.fn, "-o", umap.fn, "-t", umap.init.fn, "-c", nc, "-d",md,"-k",nn,"-n", niter)
    }
    else{
      cmd = paste("/allen/programs/celltypes/workgroups/mct-t200/Cindy_analysis/my_Tools/run_umap.py -i", rd.fn, "-o", umap.fn, "-c", nc,"-d",md,"-k",nn,"-n", niter)
    }
    print(cmd)
    system(cmd)
    alpha=1
    if(length(tmp.sampled.cells) > 50000){
      alpha=0.6
    }
    if(!is.null(meta.fields)){
      if(nc==2){
        umap.result = plot_2d_umap_anno(umap.fn, anno.df, dest.d,save.format=c("both"), meta.fields=meta.fields, show.label=TRUE,alpha=alpha)
      }
      else{
        umap.df = fread(umap.fn)
        colnames(umap.df) = c("sample_id",paste0("Dim", 1:nc))
        df = umap.df %>% left_join(anno.df)
        umap.result = plot_3d_label_multiple(umap.df, cols=meta.fields, dir=dest.d, show.label=TRUE,fn=prefix)
      }
    }                                        #plot_umap_constellation(umap.result$umap.2d, cl[sampled.cells], cl.df, select.knn.cl.df, dest.d=dest.d)    
}





find_group_umap <- function(select.cl, d, prefix, cl.sample.size=500, dims = c(2,3),nn=25, md=0.3, plot=FALSE)
{
  print(prefix)
  tmp.cl = cl[cl %in% select.cl]
  if(is.factor(tmp.cl)) {
    tmp.cl = droplevels(tmp.cl)
  }
  cl.platform = setNames(paste(tmp.cl,anno.df[names(tmp.cl),"platform"]), names(tmp.cl))
  rd.fn  = file.path(d, paste0(prefix, ".rd.parquet"))
  if(!file.exists(rd.fn)){
    fn = file.path(d,paste0(prefix,".markers.rda"))
    if(!file.exists(fn)){
      print("find markers")
      group.markers = select_markers_ds(ds, cl.bin, select.cl=select.cl,mc.cores=10)
      group.markers = intersect(group.markers, select.markers)
      save(group.markers, file=fn)
    }
    else{
      load(fn)
    }
    tmp.cells= sample_cells(cl.platform, 50)
    if(length(tmp.cells)>250000){
      tmp.cells = sample(tmp.cells, 250000)
    }    
    print("PCA")
    dat = get_cols(impute.dat.big, cols=tmp.cells, rows=group.markers)
    rd.result= rd_PCA_big(impute.dat.big, dat, select.cells=names(tmp.cl), max.dim=100, th=0, mc.cores=10)  
    tmp.cells=intersect(names(tmp.cl),comb.dat$dat.list[[1]]$col_id)  
    tmp= filter_RD(rd.result$rd.dat[tmp.cells,], rm.eigen, 0.7)
    rd.dat = rd.result$rd.dat[,colnames(tmp)]
    rd.dat.df = data.frame(sample_id=row.names(rd.dat), as.data.frame(rd.dat))
    rm(dat) 
    write_parquet(rd.dat.df, rd.fn)
    rm(rd.dat.df)
    gc()
  }
  rd.dat.df = read_parquet(file.path(d, paste0(prefix, ".rd.parquet")))
  rd.dat = as.matrix(rd.dat.df[,-1])
  row.names(rd.dat) = rd.dat.df[[1]]
  sampled.cells = sample_cells(cl.platform, cl.sample.size)
  print(length(sampled.cells))
  for(nc in dims){
    print(paste0(nc,"d.umap"))
    tmp = find_umap(unique(tmp.cl),dest.d = d, prefix=prefix, cl=tmp.cl, rd.dat=rd.dat, sampled.cells=sampled.cells, nc=nc,do.init=TRUE, nn=nn, md=md)
    umap.fn = paste0(prefix, ".umap.",nc,"d.sampled.csv")
    umap.df = as.data.frame(fread(umap.fn,header=TRUE))
    row.names(umap.df) = umap.df[[1]] 
    umap.df = umap.df[,-1]
    if(nc==2 & plot == TRUE){
      g= plot_RD_cl(umap.df, cl=cl.subclass[row.names(umap.df)], alpha=0.5, cex=0.07,fn.size=1,label.center=TRUE)
      ggsave(g, file=paste0(prefix, ".subclass.sampled.png"),height=20,width=15)
    }
    if(length(tmp.cl)> length(sampled.cells)){
      knn.result = scrattch.bigcat::get_knn_batch(dat=rd.dat[names(tmp.cl), ], ref.dat = rd.dat[sampled.cells,], k=10, method = "Annoy.Cosine", batch.size = 50000, mc.cores=20, transposed=FALSE, clear.index=TRUE, return.distance=TRUE)
      ref.cells=sampled.cells
      knn.index=knn.result$knn.index
      knn.index=knn.index[!row.names(knn.index) %in% ref.cells,]  
      w= knn.result$knn.distance
      w = w - rowMaxs(w)
      w = -w/rowMaxs(abs(w))
      w = w/rowSums(w)
      umap.other = impute_knn(knn.index, ref.cells, dat=as.matrix(umap.df), transpose_input=TRUE, transpose_output=TRUE,w=w[row.names(knn.index),])
      umap.whole = rbind(umap.df, umap.other)
    }
    else{
      umap.whole= umap.df
    }
    write.csv(umap.whole, file=file.path(d,paste0(prefix, ".umap.",nc,"d.csv")))
  }
}


plot_umap <- function(prefix,d, cex=0.3, alpha=0.7)
{
  nc=2
  print(prefix)
  umap.fn = file.path(d,paste0(prefix, ".umap.",nc,"d.csv"))
  umap.df = as.data.frame(fread(umap.fn,header=TRUE))
  row.names(umap.df) = umap.df[[1]] 
  umap.df = umap.df[,-1]
  
  if(nrow(umap.df)>300000){
    cex=0.1
    alpha=0.5
  }
  
  # cluster
  g= plot_RD_cl(umap.df, cl=cluster_id[row.names(umap.df)], alpha=alpha, cex=cex,fn.size=3,label.center=TRUE, no.shape=TRUE,
                raster=T, dpi=600)  
  #ggsave(g, file=paste0(prefix, ".cl.png"),height=20,width=15)
  ggsave(g, file=paste0(prefix, ".cl.pdf"),height=10,width=10)
  
  # supertype
  g= plot_RD_cl(umap.df, cl=cl.st[row.names(umap.df)],cl.color=st.color, alpha=alpha, cex=cex,fn.size=3,label.center=TRUE, no.shape=TRUE,
                raster=T, dpi=600)  
  #ggsave(g, file=paste0(prefix, ".cl.subclass.png"),height=20,width=15)
  ggsave(g, file=paste0(prefix, ".supertype.pdf"),height=10,width=10)
  
  # subclass == Level2
  g= plot_RD_cl(umap.df, cl=cl.level2[row.names(umap.df)],cl.color=l2.color, alpha=alpha, cex=cex,fn.size=3,label.center=TRUE, no.shape=TRUE,
                raster=T, dpi=600)  
  #ggsave(g, file=paste0(prefix, ".cl.subclass.png"),height=20,width=15)
  ggsave(g, file=paste0(prefix, ".subclass.pdf"),height=10,width=10)
  
  # class == Level1
    g= plot_RD_cl(umap.df, cl=class.cl[row.names(umap.df)],cl.color=class.col, alpha=alpha, cex=cex,fn.size=3,label.center=TRUE, no.shape=TRUE,
                raster=T, dpi=600)  
  #ggsave(g, file=paste0(prefix, ".cl.class.png"),height=20,width=15)
  ggsave(g, file=paste0(prefix, ".cl.class.pdf"),height=10,width=10)
  
  # region
  g= plot_RD_meta(umap.df, meta=factor(anno.df[row.names(umap.df),"region_label"]), meta.col=region.col, alpha=alpha, cex=cex,raster=T, dpi=600) 
  g = g + theme_void()
  #ggsave(g, file=paste0(prefix, ".region.png"),height=20,width=15)
  ggsave(g, file=paste0(prefix, ".region.pdf"),height=8,width=8)
  
  # broad region
  g= plot_RD_meta(umap.df, meta=factor(anno.df[row.names(umap.df),"broad_region"]), meta.col=bregion.col,alpha=alpha, cex=cex,legend.size=2, raster=T, dpi=600) 
  g = g + theme_void()
  #ggsave(g, file=paste0(prefix, ".broadregion.png"),height=20,width=15)
  ggsave(g, file=paste0(prefix, ".broadregion.pdf"),height=8,width=8)
  
}



