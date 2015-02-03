
#This clustering removes contigs that have reads in a representative contig!
clusterMultiMappingReads_stringent = function( contigForCountingGR_unclustered, allReads ){
  overlapMap = findOverlaps( contigForCountingGR_unclustered, allReads )  
  if( length(unique(queryHits(overlapMap))) != length(contigForCountingGR_unclustered) ){
    stop("The Clusters MUST Map all the reads provided! (This is achieved by intersection of all files for contig assembly!)as")
  }
  allClustered = TRUE
  clusteredContigs = GRanges()
  while(allClustered){
    
    if( length( overlapMap) == 0 ){
      break;
    }
    currHits = subjectHits( overlapMap[ 1 ] )
    readNamesOfContig = unique(allReads[currHits]$name)
    allReadPositions = which( allReads$name %in% readNamesOfContig )
    
    contigClusterIdx = queryHits(overlapMap)[subjectHits(overlapMap) %in% allReadPositions]
    readsToDelete = subjectHits(overlapMap)[ which(queryHits(overlapMap) %in% contigClusterIdx) ] #the other reads (with other names) that overlap this locus have to be deleted too...
    contigCluster = contigForCountingGR_unclustered[contigClusterIdx]
    toBeClusteredTo = which(contigCluster$readCount == max(contigCluster$readCount))[1]#else taking the first one currently
    
    clusteredContigs = append(clusteredContigs, contigForCountingGR_unclustered[contigClusterIdx[toBeClusteredTo]])
    
    contigForCountingGR_unclustered = contigForCountingGR_unclustered[ -contigClusterIdx  ] # delete all other locations with same read name!
    
    if(length(contigForCountingGR_unclustered) == 0 ){
      allClustered = FALSE
    } else{
      allReadsToDelete = allReads[readsToDelete]
      readToDeleteIndex = which(allReads$name %in% unique(allReads[readsToDelete]$name)) #checking for other positions of these reads from that cluster!
      allReadsToDelete = allReads[readToDeleteIndex]
      unclusteredContigsToDelete = queryHits(findOverlaps( contigForCountingGR_unclustered, allReadsToDelete ))
      if(length(unclusteredContigsToDelete) != 0){#could be that we have already captured all reads -> This is the case when all reads overlapping a cluster map to the same positions (example miRNA that maps to 3 positions)
        contigForCountingGR_unclustered = contigForCountingGR_unclustered[-unclusteredContigsToDelete]
      }
      allReads = allReads[-readToDeleteIndex]
      overlapMap = findOverlaps( contigForCountingGR_unclustered, allReads )
    }
  }
  return(clusteredContigs)
}

#
annotDFCreation = function( annotationTableSub, gffDetailColumn=9 ){
  require(R.oo)
  splittedAnnot = R.oo::trim(unlist( strsplit( annotationTableSub[,9] , ";") ))
  resVect = unlist( lapply( splittedAnnot, function(x){
    resVecttmp = unlist( strsplit( x, " " ) )
    resVect = resVecttmp[2]
    names(resVect) = resVecttmp[1]
    return(resVect)
  }) )
  
  resVect = resVect[ !duplicated(names( resVect )) ]
  
  eval( parse(text=paste0(names( resVect ), " = c()" )) )
  rm(splittedAnnot)
  
  counter = c(0)
  for( annotRow in annotationTableSub$V9 ){
    splittedAnnot = R.oo::trim(unlist( strsplit( annotRow, ";") ))
    
    resVectCurr = unlist( lapply( splittedAnnot, function(x){
      resVecttmp = unlist( strsplit( x, " " ) )
      resVect = resVecttmp[2]
      names(resVect) = resVecttmp[1]
      return(resVect)
    }) )
    
    resVectEval = resVect
    
    if( length(resVectCurr) != 0 ){
      for( i in 1:length(resVectCurr) ){
        resVectEval[names(resVectCurr)[i]] = resVectCurr[i]
      } 
      resVectEval[ !names(resVectEval) %in% names(resVectCurr)] = NA
    } else{
      resVectEval = rep(NA, length(resVect))
      names(resVectEval) = names(resVect)
    }
    
    eval( parse(text=paste0(names( resVectEval ), " = c(", names(resVectEval),",\"", resVectEval,"\")" ) ) )
    
  }
  
  eval( parse(text= paste0( "resDF = data.frame( ",paste0("\"", names(resVect), "\"=",eval(names(resVect)),collapse="," ), ")") ))
  return(resDF)
}
