###
library(igraph)
library(pwr)

###Post processing function selects most interesting links from raw ppcor results 
###according to link weight (correlation coefficient)

ppcor.post= function(PPCOR.res, coefficient.threshold){
  est= PPCOR.res$estimate
  G=graph_from_adjacency_matrix(est, mode= 'undirected', weighted = TRUE, diag= FALSE, add.colnames = NULL)
  E=get.data.frame(G, what = 'edges')
  max.value= max(E$weight, na.rm = T)
  min.value= min(E$weight, na.rm = T)
  if(max.value <=1 & min.value >=-1){valid=T}else{valid=F}
  if(valid){
    ppcor=E[which( abs(E$weight) > coefficient.threshold),]
  }else{
    ppcor=NULL
    message('PPCOR Correlation coefficients are not valid, empty network generated')
  }
  return(ppcor)
}



#Custom Rcistarget function to select high confidence links
Custom.Rcis=function(input.dir='scNET/Results',
                     pattern,
                     chosenDb,
                     output.dir='scNET/Results',
                     MinGenesetSize=0, 
                     directed=T){
  
  library(tools)
  library(RcisTarget)
  Db.dir= file_path_as_absolute(dirname(chosenDb))

  setwd(Db.dir)
   
    
  motifRankings <- importRankings(chosenDb)
  data(motifAnnotations_hgnc)
  setwd(input.dir)
  
  # read tsv files excluding other file types
  input.files=list.files(pattern = glob2rx(paste(pattern, "*.tsv", sep="")))

  Networks= lapply(input.files, read.table, header=T)
  names(Networks) <- file_path_sans_ext(basename(input.files))
  Networks= lapply(Networks, function(x){
      colnames(x)=c('gene1', 'gene2', 'weight') 
      x})
    
  TFtoTargetNetworks= lapply(Networks, function(x){
    if(directed==F){
      LinkstoKeep=c()
      LinkstoKeep.GeneStatus= data.frame(gene1= c(), gene2=c())
      for (i in 1:nrow(x)){
        Genes= as.character(x[i,1:2])
        GeneStatus= c()
        for (j in 1:length(Genes)){
          if (Genes[j] %in% colnames(getRanking(motifRankings))){
            GeneStatus[j]= 'target'
          }else{
            if (Genes[j] %in% motifAnnotations_hgnc$TF){
              GeneStatus[j]= 'TF'
            }else{GeneStatus[j]= NA}
          }
        }
        if (any(is.na(GeneStatus)) | all(GeneStatus=='target')){LinkstoKeep= LinkstoKeep}
        else{
          LinkstoKeep= c(LinkstoKeep, i)
          LinkstoKeep.GeneStatus[match(i, LinkstoKeep),1:2]= GeneStatus
        }
      }
      DbFilteredNetwork= x[LinkstoKeep,]
      for(i in 1:nrow(DbFilteredNetwork)){
        if (LinkstoKeep.GeneStatus[i,1] == 'target'){
          DbFilteredNetwork[i, 1:2]= rev( DbFilteredNetwork[i, 1:2])
        }
      }
      return(DbFilteredNetwork)
    }else{return(x)}
  })
  
  FilteredTFtoTargetNetworks= lapply(TFtoTargetNetworks, function(x){
    if(directed){
      networkFilteredforDbgenes= x[which(x$gene2 %in% colnames(getRanking(motifRankings))),]
      networkFilteredforDbgenes= networkFilteredforDbgenes[which(networkFilteredforDbgenes$gene1 %in% motifAnnotations_hgnc$TF),]
      return(networkFilteredforDbgenes)
    }else{return(x)}
  })
  
  NetworkGenesets= lapply(FilteredTFtoTargetNetworks, function(x){
    Geneset= list()
    for (i in unique(x$gene1)){
      table= x[which(x$gene1==i),]
      if (nrow(table) >= MinGenesetSize){
        Geneset[[i]]= as.character(table$gene2)
      }
    }
    return(Geneset)
  })
  
  AUC= lapply(NetworkGenesets, calcAUC, rankings= motifRankings)
  motifEnrichment= lapply(AUC, function(aucOutput){
    # Extract the TF of the gene-set name (i.e. MITF_w001):
    tf <- sapply(setNames(strsplit(rownames(aucOutput), "_"), rownames(aucOutput)), function(x) x[[1]])
    # Calculate NES and add motif annotation (provide tf in 'highlightTFs'):
    addMotifAnnotation(aucOutput, 
                       nesThreshold=3, digits=3, 
                       motifAnnot=motifAnnotations_hgnc,
                       motifAnnot_highConfCat=c("directAnnotation", "inferredBy_Orthology"),
                       motifAnnot_lowConfCat=c("inferredBy_MotifSimilarity",
                                               "inferredBy_MotifSimilarity_n_Orthology"), 
                       highlightTFs=tf)
  })
  
  motifEnrichment_selfMotifs = lapply(motifEnrichment, function(x){
    selfMotif= x[which(x$TFinDB != ""),]
    return(selfMotif)
  })
  
  EnrichmentTable= list()
  for(i in 1:length(motifEnrichment_selfMotifs)){
    EnrichmentTable[[i]]= addSignificantGenes(motifEnrichment_selfMotifs[[i]],
                                              geneSets = NetworkGenesets[[i]],
                                              rankings = motifRankings)
  }
  names(EnrichmentTable)= names(motifEnrichment_selfMotifs)
  
  
  EnrichmentTableHighConf= lapply(EnrichmentTable, function(x){
    if(nrow(x) !=0){
      HighConf= x[which(x$TFinDB=='**'),]
      return(HighConf)
    }else{return(x)}
  })
  
  HighConfRcisTargetNetwork= lapply(EnrichmentTableHighConf, function(x){
    regulons= data.frame(TF= c(), target=c())
    if(nrow(x) != 0){
      for(i in unique(x$geneSet)){
        str= x[ which(x$geneSet==i),'enrichedGenes']
        str= unlist(strsplit(str$enrichedGenes, ';'))
        str= unique(str)
        table=data.frame(TF= rep(i, length(str)), target= str)
        regulons= rbind(regulons, table)
      }
    }
    return(regulons)
  })
  
  CustomRcisTargetResult= list()
  for(i in 1:length(HighConfRcisTargetNetwork)){
    CustomRcisTargetResult[[i]]= list(AnalyzedLinks= FilteredTFtoTargetNetworks[[i]],
                                      FilteredLinks= HighConfRcisTargetNetwork[[i]])
  }
  names(CustomRcisTargetResult)= names(HighConfRcisTargetNetwork)
  
  Output.dir= file_path_as_absolute(output.dir)
  setwd(Output.dir)
  saveRDS(CustomRcisTargetResult,file = paste('RcisTarget_',pattern,'.Rds', sep = '' ))
  return(CustomRcisTargetResult)
}


#Reproducibility mesures
Intersection.index= function(netA, netB, Directed=T){
  if(nrow(netA)> nrow(netB)){
    trans=netA
    netA=netB
    netB=trans
  }
  charA= paste(netA[,1], ':', netA[,2])
  charB= paste(netB[,1], ':', netB[,2])
  if (Directed){
    inter= intersect(charA, charB)
  }else{
    charA2= paste(netA[,2], ':', netA[,1])
    charA2= c(charA, charA2)
    inter= intersect(charA2, charB)
  }
  index= length(inter)/min(length(charA), length(charB))
  return(index)
}


#Weighted Jaccard index 
wjs= function(netA, netB){
  colnames(netA)= c('gene1', 'gene2', 'weight')
  colnames(netB)= c('gene1', 'gene2', 'weight')
  #create weighted adjacency matrixes if input is an edgelist
  AG= graph.data.frame(netA, directed = T)
  Aadj= as_adjacency_matrix(AG, attr = 'weight')
  BG= graph.data.frame(netB, directed = T)
  Badj= as_adjacency_matrix(BG, attr = 'weight')
  
  #match rows and columns for both adjacency lists
  comm= intersect(row.names(Aadj), row.names(Badj))
  Badj=Badj[comm,comm]
  Aadj=Aadj[comm,comm]
  
  #ensure that everything matches
  if(all(rownames(Aadj)==rownames(Badj)) &
     all(colnames(Aadj)==colnames(Badj)) ){
    Aadj= abs(Aadj)
    Badj= abs(Badj)
    mn= pmin(Aadj, Badj)
    mx= pmax(Aadj, Badj)
    
    #calculate
    mn= sum(mn)
    mx= sum(mx)
    wjs= mn/mx
  }
  return(wjs)
}

Rcis.percent= function(RcisResult){
  values=c()
  for(i in length(RcisResult)){
    analyzed= nrow(RcisResult[[i]]$AnalyzedLinks)
    filtered= nrow(RcisResult[[i]]$FilteredLinks)
    percent= filtered/analyzed
    values=c(values, percent)
  }
  result= mean(values, na.rm=T)
  return(result)
}

Reproducibility.stats= function(Algorithm_names, Results.dir){
  Algorithm_paths <- paste(Algorithm_names, "*.tsv", sep="")
  
  Dir= file_path_as_absolute(Results.dir)
  setwd(Dir)
  res.df= data.frame(Intersection_index=rep("", length(Algorithm_names)), 
                   WJS=rep("", length(Algorithm_names)), 
                   Algorithm=rep("", length(Algorithm_names)))

  for(i in 1:length(Algorithm_names)){
    files=list.files(path= Dir, pattern= glob2rx(Algorithm_paths[i]))
    
    if (startsWith(Algorithm_names[i], "PPCOR")){
      networks= lapply(files, read.table, header=T)
    } else {networks= lapply(files, read.table)}
    index= Intersection.index(networks[[1]], networks[[2]])
    WJS= wjs(networks[[1]], networks[[2]])
    res= c(index, WJS, basename(Algorithm_names[i]))
    res.df[i,]=res
    print(paste(Algorithm_names[i], "DONE!", collapse = " : "))
     }
  return(res.df)
}


#Quantile splitting functions
library(rlist)
quantile.cut= function(net, percent){
  colnames(net)= c('gene1', 'gene2', 'weight')
  limit=quantile(net$weight, percent)
  thresh.list= list()
  for(i in 1:length(limit)){
    table= net[which(net$weight>= limit[i]),]
    thresh.list[[i]]= table
  }
  names(thresh.list)= names(limit)
  return(thresh.list)
}

quantile.stats= function(network1, network2, percent, Directed=T, label){
  FirstThresholded=  quantile.cut(network1, percent)
  SecondThresholded= quantile.cut(network2, percent)
  res.df= data.frame(network_size_net1=rep("", length(percent)), 
                     network_size_net2=rep("", length(percent)), 
                     intersection=rep("", length(percent)),
                     WJS=rep("", length(percent)), 
                     Quantile=rep("", length(percent)), 
                     Algorithm=rep("", length(percent)))
  for(i in 1:length(percent)){
    N_network1= nrow(FirstThresholded[[i]])
    N_network2= nrow(SecondThresholded[[i]])
    inter=Intersection.index(FirstThresholded[[i]], SecondThresholded[[i]], Directed= Directed)
    wjs=wjs.edge(FirstThresholded[[i]], SecondThresholded[[i]])
    res= c(N_network1, N_network2, inter, wjs, percent[i], label)
    res.df[i,]= res
  }
  return(res.df)
}

