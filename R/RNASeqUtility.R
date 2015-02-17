#' @title Summarize contig Annotation
#'
#' @description This method summarizes the contig annotation table by first checking for range based annotation then for rfam annotation then for repeat annotation 
#' @param data.frame contig annotation table (contigAnnotTable_fin.csv)
#' @return summarized contig annotation
#' @docType methods
#' @export
summarizeContigAnnotation = function(contigDF){
  currentlySupportedCols = c("contigID","g_gene_biotype","g_gene_type_ext","g_gene_name","r_repClass","r_repName","inf_targetName","inf_descriptionOfTarget")
  if( !(colnames( contigDF ) %in% currentlySupportedCols) ){
    stop( paste0("Please provide following columns in your input:\n ", paste0(currentlySupportedCols, collapse=" ") ) )    
  } 
  
  gene_type =  contigDF[,c("contigID","g_gene_biotype","g_gene_type_ext","g_gene_name")]
  gene_type$g_gene_biotype = with(gene_type, ifelse(is.na(g_gene_biotype),"none", g_gene_biotype) )
  gene_type$g_gene_name = with(gene_type, ifelse(is.na(g_gene_name),"not_annotated", g_gene_name) )
  
  #get the best feature annotation
  #start with range based annotation, then go to repeat masker if infernal is ? -> else infernal first then repeat masker
  featureAnnotType = c()
  featureAnnotName = c()
  
  featureAnnotType = with(contigDF,  ifelse( 
    #outer if statement 
    contigDF$f_featureAnnot_short == "unknown", ifelse(contigDF$inf_inc == "!" ,
                                                       sub( ".*small nucleolar RNA.*","snoRNA", sub(".*microRNA.*","miRNA",contigDF$inf_descriptionOfTarget, ignore.case=TRUE), ignore.case=TRUE ), 
                                                       ifelse(!is.na(contigDF$r_repClass), 
                                                              contigDF$r_repClass, "unknown")),
    #outer else statement
    contigDF$f_featureAnnot_short
  ) )
  
  featureAnnotName = with(contigDF,  ifelse( 
    #outer if statement 
    contigDF$f_featureAnnot_short == "unknown", ifelse(contigDF$inf_inc == "!" ,
                                                       contigDF$inf_targetName, 
                                                       ifelse(!is.na(contigDF$r_repClass), 
                                                              contigDF$r_repName, "unknown")),
    #outer else statement
    contigDF$f_feature_name
  ) )
  
  # Long form for readability but this would need to be iterated with for loop -> slow in R!:  
  #   if( contigDF$f_featureAnnot_short == "unknown" ){
  #     
  #     if( contigDF$inf_inc == "!" ){
  #       featureAnnotType = c(featureAnnotType,sub( ".*small nucleolar RNA.* ","snoRNA", sub(".*microRNA.*","miRNA",contigDF$inf_descriptionOfTarget) ))
  #       featureAnnotName = c(featureAnnotName,contigDF$inf_targetName)
  #     } else if( !is.na(contigDF$r_repClass)  ){
  #       featureAnnotType = c(featureAnnotType,contigDF$r_repClass)
  #       featureAnnotName = c(featureAnnotName,contigDF$r_repName)
  #     } else{
  #       featureAnnotType = c(featureAnnotType,"unknown")
  #       featureAnnotName = c(featureAnnotName,"unknown")
  #     }
  #     
  #   } else{
  #     featureAnnotType = c(featureAnnotType,contigDF$f_featureAnnot_short)
  #     featureAnnotName = c(featureAnnotName,contigDF$f_feature_name) 
  #   }    
  gene_type$featureAnnotType = featureAnnotType
  gene_type$featureAnnotName = featureAnnotName
  return(gene_type)
}


#' @title Randomly extract reads from fastq
#'
#' @description Extract a random subset of reads from fastq files. The files have to be in the 
#' same directory, since the method gets all fastq files from a directory and iterates through them. 
#' @param numberOfReads The number of reads one wants to extract (default: 50000)
#' @param rawDataOldDir The path to the fastq files from which the random subset is extracted
#' @param rawDataNewDir The path to the new sampled fastq files
#' @return returns TRUE
#' @docType methods
#' @export
extractRandomReadSet = function(numberOfReads=50000, rawDataOldDir, rawDataNewDir){
  require(ShortRead)
  setwd(rawDataOldDir)
  allFastq = list.files(pattern=".fastq$")
  
  tmp = lapply(allFastq, function(x){
    message("... Sampling reads ....")
    sampler <- FastqSampler(x, numberOfReads)
    set.seed(123); currSample = yield(sampler)
    currSample
    writeFastq(currSample, file.path(rawDataNewDir,x) )
    message( paste0("writing to:\n", file.path(rawDataNewDir,x) ))
  })
  return(TRUE)
}


#' @title Helper method for clustering in MappAssClust_ncRNA_Template.R object
#'
#' @description Not intended for general use, just a helper method for MappAssClust_ncRNA_Template.R
#' The method will merge the contig file generated from multiIntersect perl script with 
#' all the contigs that have been intersected with -wb option to the multiIntersected file and merged afterwards. 
#' Please see MappAssClust_ncRNA_Template.R, default inputs: multiIntersectBed_perl_CLI_cmdRes, mergeBedFile_CLI_cmdResL
#'
#' @param contigFile CmdResult object is needed as input
#' @param sampleBedFiles CmdResult object is needed as input
#' @return data.frame (Unclustered)
#' @docType methods
#' @export
##Code for Count table generation for clustering
generateCountUnclusteredTable = function( contigFile, sampleBedFiles ){
  overallContigPath = file.path( getOutFilePath(getCLIApplication(contigFile)), getOutResultName(getOutResultReference(contigFile)) )
  overallContig = read.table(overallContigPath, sep="\t", header=FALSE)
  colnames(overallContig)  = c("chr","start","end","ID","score","strand")
  overallContig = overallContig[order(overallContig$ID),]
  overallContig = overallContig[,c("chr","start","end","ID","strand")]
  
  mergedDF = do.call( cbind, mapply( function(x, y){
    mergedFilePath = file.path( getOutFilePath(getCLIApplication(x)), getOutResultName(getOutResultReference(x)) )
    
    mergedFile = read.table(mergedFilePath, sep="\t", header=FALSE)
    colnames(mergedFile)  = c("chr","start","end","ReadCount","Uniqueness","strand","ID")
    mergedFile = mergedFile[,c("ReadCount","Uniqueness","ID")]
    #Could be that there are the same IDs that have not been merged, due to the -d restriction (but the contigs are longer due to group intersection )..
    #Therefore these have to be merged -> this unfortunately may take some time since R is not the fastest with iterating
    idL = split(mergedFile, mergedFile$ID )
    mergedFile = do.call(rbind,lapply( idL, function(x){ 
      x[1,"ReadCount"] = sum( x[,1] );
      x[1,"Uniqueness"] = mean( x[,2] );
      return(x[1,]) 
      } ) )
    
    mergedTable = merge(x = overallContig, y = mergedFile[,c("ReadCount","Uniqueness","ID")], by = "ID", all.x = TRUE)
    mergedTable = mergedTable[order(mergedTable$ID),]
    
    mergedTable = mergedTable[,c("ReadCount","Uniqueness")]
    colnames( mergedTable ) = paste0(colnames( mergedTable ),"_",y)
    
    return(mergedTable)
  }, sampleBedFiles, samplesInfo$sampleName, SIMPLIFY = FALSE ))
  
  overallContig = cbind(overallContig, mergedDF)
  return(overallContig)
}

#' @title Generate ncRNA database for rnastarMapping
#'
#' @description This method downloads ensembl ncRNA data if provided with ftp adress, modifies the fasta file, 
#' subsets the fasta file to sequences of a given and generates a STAR index file if wanted. 
#'
#' @param genomeDBdir Location to where the database should be saved
#' @param ncRNAsize maximum sequence length of the fasta file (default 400nt)
#' @param nrThreads Number of threads for the rnastar program (default 6)
#' @param ensemblFTP ftp adress e.g. ftp://ftp.ensembl.org/pub/release-78/fasta/mus_musculus/ncrna/, however for convenience for mouse and human the ftp sites are stored and
#' either 'mouse' or 'human' may be used for parametrization
#' @param generateIdx Specification if an rnastar index should be created (default FALSE)
#' @return all commands generated
#' @docType methods
#' @export
#CODE FOR ncRNA DATABASE CREATION FOR RNASTAR mapping
#ftp://ftp.ensembl.org/pub/release-78/fasta/homo_sapiens/ncrna/
createSizeSelectedEnsemblncRNADB = function(genomeDBdir="/tmp/mmu", ncRNAsize=400, nrThreads=6, ensemblFTP="mouse", generateIdx=FALSE){
  dir.create(genomeDBdir)
  if( ensemblFTP == "mouse" ){
    ensemblFTP = "ftp://ftp.ensembl.org/pub/release-78/fasta/mus_musculus/ncrna/"
  }else if (ensemblFTP == "human"){
    ensemblFTP = "ftp://ftp.ensembl.org/pub/release-78/fasta/homo_sapiens/ncrna/"
  }
  cmd1 = paste0("wget --directory-prefix ",genomeDBdir," ",ensemblFTP,"*")
  message(paste0("Feching files from ensembl ", cmd1))
  system(cmd1)
  archiveFN = list.files(path=genomeDBdir, pattern = ".*.gz")
  cmd2 = paste0("gunzip ", file.path(genomeDBdir,archiveFN))
  message(paste0("Unzip file:\n", cmd2))
  system(cmd2)
  fastaFN = list.files(path=genomeDBdir, pattern = ".*.fa")
  fastaModFN = paste0(sub(".fa$","",fastaFN),".mod.fa")
  cmd3 = paste0( "python ",file.path( system.file( package="RNASeqUtility"), "data", "EliminateNewLinesFasta.py"), " -i ",file.path(genomeDBdir,fastaFN),
                 " -o ",file.path(genomeDBdir,fastaModFN) )
  message(paste0("Modifying fasta file to single line only:\n", cmd3))
  system(cmd3)
  fastaModFN_restricted = paste0(sub(".fa$","",fastaModFN),".",ncRNAsize,"nt",".fa")
  restrictedDir = file.path(genomeDBdir, paste0("restr",ncRNAsize,"nt"))
  dir.create(restrictedDir)
  cmd4 = paste0(  "awk \'!/^>/ { next } { getline seq } length(seq) <= ",ncRNAsize," { print $0 \"\\n\" seq }\' ", file.path(genomeDBdir, fastaModFN), " > ", file.path(restrictedDir, fastaModFN_restricted) )
  message(paste0("Restricting fasta to:\n",ncRNAsize," sequences ", cmd4))
  system(cmd4)
  setwd(restrictedDir)

  cmd5 = paste0( "STAR --runThreadN ", nrThreads," --runMode genomeGenerate --genomeDir ",restrictedDir," --genomeFastaFiles ",file.path(restrictedDir, fastaModFN_restricted))
  if(generateIdx){
    message(paste0("Creating IndexDB with star:\n", cmd5))
    system(cmd5)    
  }
  commands = c(cmd1,cmd2,cmd3,cmd4,cmd5)
  #return the new index directory
  return( commands )
  
}



#' @title Clustering of contigs 
#'
#' @description Contigs that contain reads that match to multiple locations are clustered to the contig with the highest read count. 
#' If more contigs have the same read count, the longest of them is chosen and then the first one in the list. 
#' Starting with the first contig (sorted by read count and length) in the unclustered list:
#' All contigs containing a read of the chosen contig are reported. 
#' If these contigs are composed by x% (readCompositionIdentity) of reads from the 
#' representative contig, they are removed from the list, including the representative contig. 
#' The representative contig is stored in the clustered list. 
#' Then the next contig of the unclustered list is chosen, until the unclustered list is empty.
#'
#' @param contigForCountingGR_unclustered GRanges object of unclustered contigs
#' @param allReads GRanges object of the single reads (obtained from each sample in bed format with bamToBed)
#' @param readCompositionIdentity The percentage of read similarity of contigs that need to be reached in order to get clustered (default 0.95)
#' When 0 is specified than one shared read leads to clustering/removing of non-representative contigs. 
#' @return clustered GRanges object
#' @docType methods
#' @export
#This clustering removes contigs that have reads in a representative contig!
clusterMultiMappingReads_stringent = function( contigForCountingGR_unclustered, allReads, readCompositionIdentity = 0.95 ){
#   tmp1 = allReads
#   tmp2 = contigForCountingGR_unclustered
#   
#   allReads = import(bamFileL[[1]], 
#                     format="BED", asRangedData=FALSE)#alignedbwt_s.bed
#   start(allReads) = start(allReads) - 1 #Strange behaviour in import of rtracklayer!
#   currBamFile = 1
  
  contigForCountingGR_unclustered = contigForCountingGR_unclustered[order(contigForCountingGR_unclustered$uniqueness, contigForCountingGR_unclustered$readCount, decreasing=TRUE)]
  overlapMap = findOverlaps( contigForCountingGR_unclustered, allReads )  
  
  if( length(unique(queryHits(overlapMap))) != length(contigForCountingGR_unclustered) ){
    warning("Some contigs have zero reads in in one sample within a group, therefore these contigs are removed from the analysis! 
            One could avoid this by intersecting all bam files with the contig files")
    contigForCountingGR_unclustered = contigForCountingGR_unclustered[ unique(queryHits(overlapMap) )]
    overlapMap = findOverlaps( contigForCountingGR_unclustered, allReads )      
  }
    
  allClustered = TRUE
  clusteredContigs = GRanges()
  while(allClustered){
    if( length( overlapMap) == 0 ){
      break;
    }
    idx = 1
    currContig = queryHits(overlapMap[idx])
    currReads = subjectHits( overlapMap[ which(queryHits(overlapMap) == currContig) ] ) #Take first entry of overlap matrix to compare
    contigReadNames = unique(allReads[currReads]$name)     #get all subject Hits (reads) for this contig
    
    allReadPositions = which( allReads$name %in% contigReadNames ) #gets the positions in the BED file (all Reads) for the reads of the contig!
    contigClusterIdx = queryHits(overlapMap)[subjectHits(overlapMap) %in% allReadPositions]#finds all other contigs containing these reads
    
    contigCluster = contigForCountingGR_unclustered[contigClusterIdx]
    toBeClusteredTo = which(contigCluster$readCount == max(contigCluster$readCount))# Fetching the contig with the highest read count!
    toBeClusteredTo = toBeClusteredTo[which( width(contigForCountingGR_unclustered[contigClusterIdx[toBeClusteredTo]]) == max( width(contigForCountingGR_unclustered[contigClusterIdx[toBeClusteredTo]]) ) )][1]#from the contigs with the highest read count taking the longest one!
    
    #Adding the contig to the set of clustered contigs!
    clusteredContigs = append(clusteredContigs, contigCluster[toBeClusteredTo])

    #contigForCountingGR_unclustered[contigClusterIdx[toBeClusteredTo]]
    
    readsFromRepContig = subjectHits(overlapMap[queryHits(overlapMap) == contigClusterIdx[toBeClusteredTo]]) #These are all reads from the representative contig!
    repContigReadNames = unique(allReads[readsFromRepContig]$name) 
    allReadPositionsRepContig = which( allReads$name %in% repContigReadNames )
    #From the representative contig all reads are obtained and eliminated from the overlapMap (therefore if another contig would consist 100% of the same reads it will get eliminated)
    #Otherwise the contig will be kept and may be of interest -> one could delete some further contigs by x percent...
    
    overlapMap_filt =  overlapMap[which( !subjectHits(overlapMap) %in% allReadPositionsRepContig)]#This would be the filtered overlapMap!
    #check for contig leftovers...
    #Filtering by 95% read composition identity!
    leftoverContigs = table( queryHits(overlapMap_filt[ queryHits( overlapMap_filt ) %in% unique(contigClusterIdx)]))
    filtered_leftoverContigs = names(leftoverContigs[ifelse( 1 - (leftoverContigs / (contigCluster[toBeClusteredTo]$readCount + leftoverContigs) ) >= readCompositionIdentity , TRUE, FALSE )]) #left with less than 95% of the reads in main cluster!
    overlapMap_filt = overlapMap_filt[which(!queryHits(overlapMap_filt) %in% filtered_leftoverContigs)]
        
    overlapMap = overlapMap_filt
  }
  return(clusteredContigs)
}








#' @title Deprecated version of clustering
#'
#' @description This method is deprecated do not use!
#' @param contigForCountingGR_unclustered
#' @param allReads
#' @return clusteredContigs
#' @docType methods
#' @export
#This clustering removes contigs that have reads in a representative contig!
clusterMultiMappingReads_stringent_deprecated = function( contigForCountingGR_unclustered, allReads ){
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
