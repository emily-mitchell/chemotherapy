### Add signatures onto trees

install.packages("ape")

exposures <- read.table("exposures_new.txt", header=T,check.names =F, sep="\t",quote = "", row.names = 1)

rownames(exposures)=exposures[,1]
exposures <- exposures[,-1]


#select colonies samples only
exposures_nanoseq = exposures[2662:2740,]

exposures_colonies = exposures[1:2661,]

exposures_colonies <- t(exposures_colonies)

patients=read.table("patients.txt")[,1]

sigs=rownames(exposures_colonies)


#sig_profiles=mut_example_multi@comp_categ_distn$mean

br_lengths <- NULL
for (patient in patients){
  if(file.exists(paste0("/trees/",patient,".tree"))){
    tree = read.tree(paste0("/trees/",patient,".tree"))
    tree_df=fortify(tree)
    
    cols= c("wheat2","powderblue","darkslateblue","maroon3","#ffac3c","darkorange4","springgreen3","#0b6623","darkolivegreen2","mediumorchid3","pink","red","ivory4")
    
    #cols= c("wheat2","powderblue","darkslateblue","maroon3","#ffac3c","darkorange4","springgreen3","#0b6623","darkolivegreen2","mediumorchid3","pink","red2", "ivory4")
    
    names(cols)=sigs
    samples=colnames(exposures_colonies)[grepl(patient,colnames(exposures_colonies))]
    branches=sub('.*_', '', samples)
    
    pdf(paste0("/trees/",patient,"_tree_signatures_small_new.pdf"), height=4)
    plot(tree,show.tip.label=F, cex= 0.4)
    for (k in 1:length(samples)){
      n=as.numeric(branches[k])
      x_end=tree_df$x[n]
      x_start=tree_df$x[tree_df$parent[n]]
      x_intv=x_end-x_start
      y=node.height(tree)[n]
      l=tree_df$branch.length[n]
      t=samples[k]
      br_lengths <- rbind(br_lengths,cbind(t, l))
      tipnum=sum(tree_df$isTip)
      title(main=patient)
      for (s in sigs){
        x_end=x_start+ exposures_colonies[s,samples[k]]*x_intv
        rect(ybottom=y-min(0.02*tipnum,0.4),ytop=y+min(0.02*tipnum,0.4),xleft=x_start,xright=x_end,col=cols[s], border = cols[s] , lwd = 0.5)
        x_start=x_end
        
      }
    }
    
    axisPhylo(side = 1,backward=F)
    #      legend("topleft",title="Signatures", legend=paste0(sigs), 
    #           fill=cols, bty="n",cex=0.4, ncol=1, xjust=0.5)
    
    #axisPhylo(side = 1,backward=F) 
    #legend(x = -1.1, y = max(y), title = "Signatures", legend = paste0(sigs), 
    #fill = cols, bty = "n", cex = 0.6, ncol = 1, xjust = 1, xpd = TRUE)
    dev.off()
    
  }
}