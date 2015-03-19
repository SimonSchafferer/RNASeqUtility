
#' @title Calculate Normalized Read Count Coverage
#'
#' @description This method returns a data frame containing normalized coverage values for each position. It requires the output of calcReadCountCoveragePerCondition or generateCoverageDF_contigs as data.frame or list
#' @param coverageL  output of calcReadCountCoveragePerCondition/generateCoverageDF_contigs
#' @param deseqds DESeq object of the current data
#' @return data.frame containing the normalized read counts for each position in order to create a histogram of read counts for each ncRNA
#' @docType methods
#' @export
calcNormReadCountCoverage = function( coverageL, deseqds ){
    #Converting list to data.frame
    coverageDF = do.call(rbind, coverageL)
  
  if( sum(colnames(coverageDF) == "Sample") == 0 ){
    stop(paste0("Column 'Sample' must be present in coverageDF!"))
  }
  
  coverageDFL = split(coverageDF, factor( coverageDF$Sample, levels=unique(coverageDF$Sample), ordered=FALSE  ) )
  if( sum(names(coverageDFL) == dds$sampleName) != length(dds$sampleName) ){
    stop("SampleNames must be the same in the DESeq object and the coverageDF")
  }
  
  coverageDF = do.call(rbind,mapply(function(x,y){
    x$ReadCountNorm = x$ReadCount / y
    return(x)
  }, coverageDFL, sizeFactors(deseqds), SIMPLIFY = FALSE) )
  
  return(coverageDF)
}

#' @title Splits coverage list per condition
#'
#' @description Splits the list obtained by calcReadCountCoveragePerCondition/generateCoverageDF_contigs by condition and calculates the mean read coverage per condition
#' @param coverageL  output of calcReadCountCoveragePerCondition/generateCoverageDF_contigs
#' @param deseqds DESeq object of the current data
#' @param samplesInfo describing the experiment (produced by configuration template)
#' @return data.frame containing the normalized read counts for each position in order to create a histogram of read counts for each ncRNA
#' @docType methods
#' @export
splitReadCountCoveragePerCondition = function( coverageL, deseqds, samplesInfo ){
  
  coverageDF = RNASeqUtility::calcNormReadCountCoverage(coverageL, deseqds)
  #Mapping condition to sample groups
  coverageDF = do.call(rbind, lapply(  split(coverageDF,coverageDF$Sample), function(x){
    x$condition = samplesInfo[ which(samplesInfo$sampleName == x$Sample[1]),]$condition
    return(x)
  }))
  
  #Splitting the data frame by condition and calculating the mean coverage over samples
  coverageDFL = split( coverageDF, coverageDF$condition)
  coverageDFL = lapply(coverageDFL, function(x){
    
    ncRNAClassL = split(x, x$ID)
    
    ncRNAClassL_mod = lapply( ncRNAClassL , function(y){
      coordL = split(y, y$Idx)
      
      coordL_modL = lapply( coordL, function(z){
        z$ReadMeans = mean( z$ReadCount, na.omit=TRUE )
        z$ReadMeansNorm = mean( z$ReadCountNorm, na.omit=TRUE )
        return(z)
      })
      coordL_mod = do.call(rbind,coordL_modL)
      
      return(coordL_mod)
    })
    return(do.call(rbind, ncRNAClassL_mod))
  })
  
  return( do.call(rbind, coverageDFL) )
  
}


#' @title Generate coverage data frame for specific candidates
#'
#' @description This method returns a data frame containing coverage values for each position. It is intended for genome mapping results. For ncRNA mappings please see: generateCoverageDF_ncRNAmappings
#' @param contigsGR GRanges file containing the coordinates of the contigs and their contigID (specified as contigID in the metadata)
#' @param contigResTable Selected candidates from the result of differential expression analysis from DESeq2 combined with the minimum annotation dataframe
#' @param bedgraphMapping_plusL the cmdResult object of the bedgraph plus mapping (genomeCovNcRNAMapping_plusL) 
#' @param bedgraphMapping_minusL the cmdResult object of the bedgraph minus mapping (genomeCovNcRNAMapping_minusL) 
#' @param sampleNames names of the samples from samplesInfo table
#' @param annotationGR An otpional annotation GRanges object, when provided the annotation will be used as reference interval
#' @return data.frame containing the read counts for each position in order to create a histogram of read counts for each ncRNA
#' @docType methods
#' @export
generateCoverageDF_contigs = function(contigsGR, contigResTable, bedgraphMapping_plusL, bedgraphMapping_minusL, sampleNames, annotationGR){
  require(rtracklayer)
  require(GenomicRanges)
  
  colsRequired = c("ID","UID","biotype","Name","log2FoldChange","pvalue","padj","strandness")
  if( sum( colsRequired %in% colnames(resAnnot) ) != length(colsRequired) ){
    stop( paste0( "Columns: ",paste0(colsRequired,collapse=", ")," have to be provided!"))
  }
  
  #Combining to GRanges object
  roisContigs = contigsGR[ contigsGR$contigID %in% contigResTable$UID ]
  elementMetadata(roisContigs) = merge( elementMetadata( roisContigs ), contigResTable, by.x="contigID", by.y="UID" ,sort=FALSE)
  roisContigsL = GenomicRanges::split(roisContigs, roisContigs$ID)
  
  bedgraphncRNAMapping_plusL = lapply( bedgraphMapping_plusL, function(x){
    rtracklayer::import(file.path(getOutFilePath(getCLIApplication(x)), getOutResultName(getOutResultReference(x))), 
           asRangedData=FALSE, format="bedGraph")
  } )
  bedgraphncRNAMapping_minusL = lapply( bedgraphMapping_minusL, function(x){
    rtracklayer::import(file.path(getOutFilePath(getCLIApplication(x)), getOutResultName(getOutResultReference(x))), 
           asRangedData=FALSE, format="bedGraph")
  } )
  
  
  if(missing(annotationGR)){
    
    roiL = vector("list", length(bedgraphncRNAMapping_plusL))
    
    for(i in c(1:length(bedgraphncRNAMapping_plusL))){
      bgminus = bedgraphncRNAMapping_minusL[[i]]
      bgplus = bedgraphncRNAMapping_plusL[[i]]
      
      roiL_contig = lapply( roisContigsL  , function(x){
        
        if( as.character(strand(x)) == "+"  ){
          rnaOI = GenomicRanges::subsetByOverlaps(bgplus,x )
        }else{
          rnaOI = GenomicRanges::subsetByOverlaps(bgminus,x )
        }    
        
        if(length(rnaOI) != 0){
          #Filling the gaps in the coverage file by creating a one base based annotation with zero score! 
          rnaOI_fill = GRanges( as.character(seqnames(rnaOI))[1], 
                                IRanges(
                                  start=start(rnaOI)[1]:(end(rnaOI)[length(rnaOI)]-1),
                                  end=(start(rnaOI)[1]+1):end(rnaOI)[length(rnaOI)] ),
                                score=rep(0, end(rnaOI)[length(rnaOI)]-start(rnaOI)[1] ) )
          
          fillGR = GenomicRanges::setdiff( rnaOI_fill, rnaOI )
          
          if( length(fillGR) > 0 ){
            fillGR$score = 0
            rnaOI = IRanges::append( rnaOI, fillGR)
          }          
        }
        
        rnaOICov = rep(rnaOI$score, width(rnaOI))#
        if( length(rnaOICov) == 0  ){
          rnaOICovDF = NULL
        } else{
          rnaOICovDF = data.frame( "Idx"=1:length(rnaOICov), 
                                   "ReadCount"=rnaOICov, 
                                   "group"=x$biotype, 
                                   "Name"=x$Name,
                                   "ID"=x$ID, 
                                   "regulation"=ifelse( x$log2FoldChange > 0, "up", "down"), 
                                   "log2FoldChange"=x$log2FoldChange,
                                   "pvalue"=x$pvalue,
                                   "padj"=x$padj,
                                   "strandness"=x$strandness
          )
          
        } 
        
        return(rnaOICovDF)
      } )
      roi_ncOtherDF = do.call(rbind, roiL_contig)
      roi_ncOtherDF$Sample = ifelse(missing(sampleNames),paste0("Sample",i), sampleNames[i])
      roiL[[i]] = roi_ncOtherDF
    }
    return(roiL)
  } else{
    
    currMap = GenomicRanges::findOverlaps(annotationGR, roisContigs)
    currMap = currMap[ !queryHits(currMap) %in% queryHits(currMap)[which(duplicated( queryHits(currMap)))]  ]
    
    roisContigs = roisContigs[subjectHits(currMap)]
    annotationGR = annotationGR[queryHits(currMap)]
    
    roiL = vector("list", length(bedgraphncRNAMapping_plusL))
    
    for(i in c(1:length(bedgraphncRNAMapping_plusL))){
      bgminus = bedgraphncRNAMapping_minusL[[i]]
      bgplus = bedgraphncRNAMapping_plusL[[i]]
      
      roicL = vector("list", length(roisContigs)) 
      for(k in 1:length(roisContigs)){
        contig = roisContigs[k]
        refAnnot = annotationGR[k]
        
        if( as.character(strand(contig)) == "+"  ){
          rnaOI = GenomicRanges::subsetByOverlaps(bgplus,refAnnot )
        }else{
          rnaOI = GenomicRanges::subsetByOverlaps(bgminus,refAnnot )
        }    
        
        if(length(rnaOI) == 0){
          #No coverage information available in this library!
          return(data.frame())
        }
        
        if( start(rnaOI)[1] > start(refAnnot) ){
          toAdd = (length( start(refAnnot):start(rnaOI)[1] ))
          beforeToAdd = GRanges( seqnames(rnaOI)[1], IRanges( (start(rnaOI)[1]-toAdd+1), (start(rnaOI)[1]-1) ), score=0 )
          rnaOI = IRanges::append(rnaOI, beforeToAdd)
        } else{
          start(rnaOI)[1] = start(refAnnot)
        }
        
        if(end(refAnnot) > end(rnaOI)[length(rnaOI)]){
          toAdd = length( end(rnaOI)[length(rnaOI)]:end(refAnnot) )
          afterToAdd = GRanges( seqnames(rnaOI)[1], IRanges( (end(rnaOI)[length(rnaOI)]+1), 
                                                             (end(rnaOI)[length(rnaOI)]+toAdd-1  ) ), score=0 )
          rnaOI = IRanges::append(rnaOI, afterToAdd)
        } else{
          end(rnaOI)[length(end(rnaOI))] = end(refAnnot)
        } 
        
        #Filling the gaps in the coverage file by creating a one base based annotation with zero score! 
        refAnnot_fill = GRanges( as.character(seqnames(refAnnot))[1], 
                                 IRanges(
                                   start=start(refAnnot):(end(refAnnot)-1),
                                   end=(start(refAnnot)+1):end(refAnnot) ),
                                 score=rep(0,width(refAnnot)-1))
        
        fillGR = GenomicRanges::setdiff( refAnnot_fill, rnaOI )
        
        if( length(fillGR) > 0 ){
          fillGR$score = 0
          rnaOI = IRanges::append( rnaOI, fillGR)
        }
        
        rnaOI = sort(rnaOI)
        
        rnaOICov = rep(rnaOI$score, width(rnaOI))#
        
        rnaOICovDF = data.frame( "Idx"=1:length(rnaOICov), 
                                 "ReadCount"=rnaOICov, 
                                 "group"=contig$biotype, 
                                 "Name"=contig$Name,
                                 "ID"=contig$ID, 
                                 "regulation"=ifelse( contig$log2FoldChange > 0, "up", "down"), 
                                 "log2FoldChange"=contig$log2FoldChange,
                                 "pvalue"=contig$pvalue,
                                 "padj"=contig$padj,
                                 "strandness"=contig$strandness
        )
        roicL[[k]] = rnaOICovDF
      }
      roiDF = do.call(rbind, roicL)
      roiDF$Sample = ifelse(missing(sampleNames),paste0("Sample",i), sampleNames[i])
      roiL[[i]] = roiDF
    }
    return(roiL)
  }
}


#' @title Generate coverage data frame for specific candidates
#'
#' @description This method returns a data frame containing coverage values for each position. It is intended for ncRNAmapping results. For contigs please see: generateCoverageDF_contigs
#' @param ncRNAMappingsDF Selected candidates from the result of differential expression analysis from DESeq2 combined with the minimum annotation dataframe
#' @param bedgraphMapping_plusL the cmdResult object of the bedgraph plus mapping (genomeCovNcRNAMapping_plusL) 
#' @param bedgraphMapping_minusL the cmdResult object of the bedgraph minus mapping (genomeCovNcRNAMapping_minusL) 
#' @param sampleNames names of the samples from samplesInfo table
#' @return data.frame containing the read counts for each position in order to create a histogram of read counts for each ncRNA
#' @docType methods
#' @export
generateCoverageDF_ncRNAmappings = function( ncRNAMappingsDF, bedgraphMapping_plusL, bedgraphMapping_minusL, sampleNames ){
  require(rtracklayer)
  require(GenomicRanges)
  colsRequired = c("ID","UID","biotype","Name","log2FoldChange","pvalue","padj","strandness")
  if( sum( colsRequired %in% colnames(resAnnot) ) != length(colsRequired) ){
    stop( paste0( "Columns: ",paste0(colsRequired,collapse=", ")," have to be provided!"))
  }
  
  ncRNAMappingsDF$ID = sub("_plus.*$","",ncRNAMappingsDF$UID)
  ncRNAMappingsDF$ID = sub("_minus.*$","",ncRNAMappingsDF$ID)
  
  bedgraphncRNAMapping_plusL = lapply( bedgraphMapping_plusL, function(x){
    rtracklayer::import(file.path(getOutFilePath(getCLIApplication(x)), getOutResultName(getOutResultReference(x))), 
           asRangedData=FALSE, format="bedGraph")
  } )
  bedgraphncRNAMapping_minusL = lapply( bedgraphMapping_minusL, function(x){
    rtracklayer::import(file.path(getOutFilePath(getCLIApplication(x)), getOutResultName(getOutResultReference(x))), 
           asRangedData=FALSE, format="bedGraph")
  } )
  
  ncRNAMappingsDFL = split(ncRNAMappingsDF, ncRNAMappingsDF$UID)
  
  roiL = vector("list", length(bedgraphncRNAMapping_plusL))
  
  for(i in c(1:length(bedgraphncRNAMapping_plusL))){
    bgminus = bedgraphncRNAMapping_minusL[[i]]
    bgplus = bedgraphncRNAMapping_plusL[[i]]
    
    roi_ncOtherL = lapply( ncRNAMappingsDFL  , function(x){
      
      #First decide if plus or minus according to UID
      if( length(grep("minus",x$UID)) == 0  ){
        currGR = bgplus[seqnames( bgplus ) ==  x$ID]
      }else{
        currGR = bgminus[seqnames( bgminus ) ==  x$ID]
      }    
      
      if(length(currGR) == 0) {
        #Too less reads, in the mapping file are no reads available
        return(data.frame())
      }
      
      #Filling the gaps in the coverage file by creating a one base based annotation with zero score! 
      rnaOI_fill = GRanges( as.character(seqnames(currGR))[1], 
                            IRanges(
                              start=start(currGR)[1]:(end(currGR)[length(currGR)]-1),
                              end=(start(currGR)[1]+1):end(currGR)[length(currGR)] ),
                            score=rep(0, end(currGR)[length(currGR)]-start(currGR)[1] ) )
      
      fillGR = GenomicRanges::setdiff( rnaOI_fill, currGR )
      
      if( length(fillGR) > 0 ){
        fillGR$score = 0
        currGR = IRanges::append( currGR, fillGR)
      }
      
      rnaOICov = rep(currGR$score, width(currGR))#
      rnaOICovDF = data.frame( "Idx"=1:length(rnaOICov), 
                               "ReadCount"=rnaOICov, 
                               "group"=x$biotype, 
                               "Name"=x$Name,
                               "ID"=x$ID, 
                               "regulation"=ifelse( x$log2FoldChange > 0, "up", "down"), 
                               "log2FoldChange"=x$log2FoldChange,
                               "pvalue"=x$pvalue,
                               "padj"=x$padj,
                               "strandness"=x$strandness )
      return(rnaOICovDF)
    } )
    roi_ncOtherDF = do.call(rbind, roi_ncOtherL)
    roi_ncOtherDF$Sample = ifelse(missing(sampleNames),paste0("Sample",i), sampleNames[i])
    roiL[[i]] = roi_ncOtherDF
  }
  return(roiL)
}


#' @title Generate Count table
#'
#' @description This method combines the read counts after counting them  with mergeBed from bam files
#' @param bamToBedAndMergencRNAL bed files after mapping
#' @return read count table
#' @docType methods
#' @export
generateCountFromMappingDF = function( bamToBedAndMergencRNAL, samplesInfo ){
  require(rtracklayer)
  require(GenomicRanges)
  
  ncRNAReadCountL = lapply( bamToBedAndMergencRNAL, function(x){
    x = x[[2]]
    dir = getOutFilePath(getCLIApplication(x))
    return( rtracklayer::import( file.path( dir,getOutResultName(getOutResultReference(x))), format="BED",asRangedData=FALSE ) )
  })
  
  #Create a reference dataframe containing all entries (Ensembl ids plus strand as ID)
  overallTable = unique(do.call( rbind, lapply( ncRNAReadCountL, function(x){
    xdf = GenomicRanges::as.data.frame(x)
    colnames(xdf) = c("transcript_id", "mapStart","mapEnd","mapWidth","mapStrand","readCount","mapQ")
    xdf$UID = paste0(xdf$transcript_id, ifelse(xdf$mapStrand == "-", "_minus", "_plus") )
    return(xdf[,c("transcript_id","UID")])
  }) ) )
  rownames(overallTable) = 1:dim(overallTable)[1]
  
  
  #Now merge with counted reads
  readsCountDF = do.call( cbind, lapply( ncRNAReadCountL, function(x){
    xdf = GenomicRanges::as.data.frame(x)
    colnames(xdf) = c("transcript_id", "mapStart","mapEnd","mapWidth","mapStrand","readCount","mapQ")
    xdf$UID = paste0(xdf$transcript_id, ifelse(xdf$mapStrand == "-", "_minus", "_plus") )
    #Merge The duplicated UIDs
    xdf = do.call(rbind, lapply( split(xdf, xdf$UID), function(x){
      xnew = x[1,]
      xsize = dim(x)[1]
      if( xsize != 1 ){
        xnew$mapStart = min(as.numeric(x$mapStart))
        xnew$mapEnd = min(as.numeric(x$mapEnd))
        xnew$readCount = sum(as.numeric(x$readCount),na.rm = TRUE)
        xnew$mapQ = mean(as.numeric(x$mapQ),na.rm = TRUE)
      }
      return(xnew)
    }))
    xdf_merged = merge(overallTable, xdf[,c("UID", "readCount","mapStart","mapEnd")], by="UID", all.x=TRUE)
    return(xdf_merged)
  } ) )
  
  #Change for mapStart take minimum and map End take maximum? -> else rowMeans should be ok too, but round!
  mapStart = readsCountDF[,grep("mapStart", colnames(readsCountDF)) ]
  mapStart = apply(mapStart, 1, min, na.rm=TRUE)
  mapEnd = readsCountDF[,grep("mapEnd", colnames(readsCountDF)) ]
  mapEnd = apply(mapEnd, 1, min, na.rm=TRUE)
  
  
  readsCountDF = readsCountDF[, c(1, grep("readCount", colnames(readsCountDF))) ]
  colnames(readsCountDF)[1] = "UID"
  colnames(readsCountDF)[2:length(colnames(readsCountDF))] = as.character(samplesInfo$sampleName)
  #Changing all NA values to zero
  readsCountDF[,2:length(colnames(readsCountDF))] = apply(readsCountDF[,2:length(colnames(readsCountDF))], 2, function(x){
    ifelse( is.na(x), 0, x)
  })
  readsCountDF$mapStart = mapStart
  readsCountDF$mapEnd = mapEnd
  
  return( readsCountDF )
  
}

#' @title Summarize contig Annotation
#'
#' @description This method summarizes the contig annotation table by first checking for range based annotation then for rfam annotation then for repeat annotation 
#' @param data.frame contig annotation table (contigAnnotTable_fin.csv)
#' @return summarized contig annotation
#' @docType methods
#' @export
summarizeContigAnnotation = function(contigDF){
  currentlySupportedCols = c("contigID","g_gene_biotype","g_gene_type_ext","g_gene_name","r_repClass","r_repName","inf_targetName","inf_descriptionOfTarget")
  if( sum(colnames( contigDF ) %in% currentlySupportedCols) == length(colnames( contigDF )) ){
    stop( paste0("Please provide following columns in your input:\n ", paste0(currentlySupportedCols, collapse=" ") ) )    
  } 
  
  gene_type =  contigDF[,c("contigID","g_gene_biotype","g_gene_type_ext","g_gene_name")]
  gene_type$g_gene_biotype = with(gene_type, ifelse(is.na(g_gene_biotype),"none", g_gene_biotype) )
  gene_type$g_gene_name = with(gene_type, ifelse(is.na(g_gene_name),"not_annotated", g_gene_name) )
  
  #get the best feature annotation
  #start with range based annotation, then go to repeat masker if infernal is ? -> else infernal first then repeat masker
  featureAnnotType = c()
  featureAnnotName = c()
  annotationType = c()
  
  featureAnnotType = with(contigDF,  ifelse( 
    #outer if statement 
    contigDF$f_featureAnnot_short == "unknown", ifelse(!is.na(contigDF$inf_inc) & contigDF$inf_inc == "!" ,
                                                       sub( ".*small nucleolar RNA.*","snoRNA", sub(".*microRNA.*","miRNA",contigDF$inf_descriptionOfTarget, ignore.case=TRUE), ignore.case=TRUE ), 
                                                       ifelse(!is.na(contigDF$r_repClass), 
                                                              contigDF$r_repClass, "unknown")),
    #outer else statement
    contigDF$f_featureAnnot_short
  ) )
  
  featureAnnotName = with(contigDF,  ifelse( 
    #outer if statement 
    contigDF$f_featureAnnot_short == "unknown", ifelse(!is.na(contigDF$inf_inc) & contigDF$inf_inc == "!" ,
                                                       contigDF$inf_targetName, 
                                                       ifelse(!is.na(contigDF$r_repClass), 
                                                              contigDF$r_repName, "unknown")),
    #outer else statement
    contigDF$f_feature_name
  ) )

  
  annotationType = with(contigDF,  ifelse( 
    #outer if statement 
    contigDF$f_featureAnnot_short == "unknown", ifelse(!is.na(contigDF$inf_inc) & contigDF$inf_inc == "!" ,
                                                       "infernal", 
                                                       ifelse(!is.na(contigDF$r_repClass), 
                                                              "repeatMaskr", "unknown")),
    #outer else statement
    "featureAnnot"
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
  gene_type$featureAnnotType = featureAnnotType#ifelse( is.na(featureAnnotType), "none",annotationType)
  gene_type$featureAnnotName = featureAnnotName
  gene_type$annotationType = annotationType#ifelse( is.na(annotationType), "none",annotationType)
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
  overlapMap = GenomicRanges::findOverlaps( contigForCountingGR_unclustered, allReads )  
  
  if( length(unique(queryHits(overlapMap))) != length(contigForCountingGR_unclustered) ){
    warning("Some contigs have zero reads in in one sample within a group, therefore these contigs are removed from the analysis! 
            One could avoid this by intersecting all bam files with the contig files")
    contigForCountingGR_unclustered = contigForCountingGR_unclustered[ unique(queryHits(overlapMap) )]
    overlapMap = GenomicRanges::findOverlaps( contigForCountingGR_unclustered, allReads )      
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
