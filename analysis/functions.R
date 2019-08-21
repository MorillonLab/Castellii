##  ***************************************************
##  *  Ving's methods                                 *
##  *  http://vm-gb.curie.fr/ving/                    *
##  *  Descrimes et al. BMC Res Notes. 2015; 8: 419.  * 
##  ***************************************************

library("multicore")

## Converts the flag numbers into binary 
integer.base.b=function(x, b=2) {
  xi=as.integer(x)
  if(any(is.na(xi) | ((x-xi)!=0)))
    print(list(ERROR="x not integer", x=x))
  N=length(x)
  xMax=max(c(x,1))
  ndigits=11
  Base.b=array(NA, dim=c(N, ndigits))
  for(i in 1:ndigits) {#i=1
    Base.b[, ndigits-i+1]=(x %% b)
    x=(x %/% b)
  }
  Base.b
} 

## Combines two data frame 
listOfDataFrame2DataFrame=function(x,vertical=TRUE) {
  if(length(x)==1) {
    x[[1]]
  }
  else {
    if(vertical) {
      ncol=dim(x[[1]])[2]
      vect=unlist(lapply(x,t))
      retour=as.data.frame(t(array(vect,dim=c(ncol,length(vect)/ncol))),stringsAsFactors=FALSE)
      names(retour)=names(x[[1]])
    }
    else {
      nline=dim(x[[1]])[1]
      vect=unlist(x)
      retour=as.data.frame(array(vect,dim=c(nline,length(vect)/nline)))
      names(retour)=unlist(lapply(x,names))
    }
    retour
  }
}

extractFromBam=function(file,which,what) {
  return(scanBam(file, param=ScanBamParam(which=which,what=what))[[1]][[1]])
}

## Returns a list from cigar expression
interpreteCIGAR=function(cigar) {
  cigar_un=strsplit(unique(cigar),split="")
  n_cigar_un=length(cigar_un)
  taille_cigar=list()
  analise_cigar=function(cigar_) {
    cigar_sortie=list()
    acc=""
    for(j in 1:length(cigar_)) {
      if(sum(cigar_[j]==as.character(0:9))==1) {
	acc=paste(acc,cigar_[j],sep="")
      }
      else {
	cigar_sortie[[length(cigar_sortie)+1]]=as.integer(acc)
	cigar_sortie[[length(cigar_sortie)+1]]=cigar_[j]
	acc=""
      }
    }
    return(cigar_sortie)
  }
  cigar_interprete=lapply(cigar_un,analise_cigar)
  names(cigar_interprete)=unique(cigar)

  return(cigar_interprete)
}

# prend un CIGAR splités et retourne la taille occupé par le read sur la séquence génomique (introns compris)
calcule_longueur_cigar=function(cigar_) {
  lon=0
  N=length(cigar_)
  for(j in seq(2,N,2)) {
    if(cigar_[[j]]!="I") {
      lon=lon+cigar_[[j-1]]
    }
  }
  return(lon)
}

# prend un CIGAR splités et retourne les positions 
calcule_junction_cigar=function(cigar_) {
  retour=list()
  lon=0
  N=length(cigar_)
  for(j in seq(2,N,2)) {
    if(cigar_[[j]]!="I") {
      lon=lon+cigar_[[j-1]]
    }
    if(cigar_[[j]]=="N") {
      retour[[length(retour)+1]]=c(lon-cigar_[[j-1]]+1,lon)
    }
  }
  return(retour)
}

## Returns a list of numbers of single read with their coordinates
compresse_coordonnees=function(debut,fin) {
  if(length(debut)==0) {
    return(list(numeric(),numeric(),numeric()))
  }
  else {
    tmp=sort(paste(debut,fin,sep="_"))
    tmp_rle=rle(tmp)
    poids=tmp_rle$lengths
    values_split=strsplit(tmp_rle$values,split="_")
    doit=function(j) {
      return(as.integer(values_split[[j]][1]))
    }
    debut_uni=sapply(1:length(poids),doit)
    doit=function(j) {
      return(as.integer(values_split[[j]][2]))
    }
    fin_uni=sapply(1:length(poids),doit)
    ordre_debut=order(debut_uni)
    return(list(debut_uni[ordre_debut],fin_uni[ordre_debut],poids[ordre_debut]))  
  }
}

RDataFileName=function(file) {
  return(paste(file,".RData",sep=""))
}

## Function converts and extracts the infos from bamfile to bamHandler
# bamHandler : data structure to describe a bam file, using 3 (or 4) vectors per chromosome
# for a given sample : 
# the first vector gives the start position (sorted) of all reads in the bam,
# the second vector gives the stop position of the corresponding reads
# the third vector gives the number of reads with coordinates corresponding to start & stop values given in first and second vectors
# the fourth vector is used for normalisation, and correspond to the third vector multiplied by normalization coefficient
# same is done for exon junction
# also contains information on chromosomes (size) and mapping stats
readBam_=function(file_,insert_max_=2000,stranded_=TRUE,ncore_=1,libraryType_=c("standard","inverse"),fileNameRData_=NA,normalized_=NULL,chrName_=NULL,from_=1,to_=NULL) {
  suppressPackageStartupMessages(require("Rsamtools"))
  suppressPackageStartupMessages(require("GenomicRanges"))
## Declaration of variables
  flagstat=numeric(11)
  names(flagstat)=c("total","duplicates","mapped","paired","read1","read2","properly paired","itself and mate mapped","singletons","mate mapped on a different chr","QC-failed")
  genome_info=scanBamHeader(file_)[[1]]$targets
  noms_chromosomes=names(genome_info)
  longueur_chromosomes=as.integer(genome_info)
  nombre_chromosomes=length(noms_chromosomes)
  brin_F=list()
  brin_R=list()
  brin_F_junction=list()
  brin_R_junction=list()
  pas=c(1,2,6,7)
  i_zone=0
  if(is.null(chrName_)) {
    chrName__=noms_chromosomes
  }else {
    chrName__=chrName_
  }
## Fragments identification
  for(i in (1:nombre_chromosomes)) {
    i_zone=i_zone +1
    nom_chromo=noms_chromosomes[i]
    lon_chromo=longueur_chromosomes[i]
    
    if(!(nom_chromo %in% chrName__)) {
      brin_F[[i]]=list(numeric(),numeric(),numeric())
      brin_R[[i]]=list(numeric(),numeric(),numeric())
      brin_F_junction[[i]]=list(numeric(),numeric(),numeric())
      brin_R_junction[[i]]=list(numeric(),numeric(),numeric())
    }else {
      if(is.null(to_)) {
	to_i=lon_chromo
      }else {
	to_i=to_[min(i_zone,length(to_))]
      }
      from_i=from_[min(i_zone,length(from_))]
      
      commande=paste("RangesList(`",nom_chromo,"`=IRanges(",from_i,",",to_i,"))",sep="")
      expr=try(parse(text=commande),TRUE)
## Function used from GenomicRanges package
      which=eval(expr)
      what=c("flag","mpos","cigar","mrnm","isize")
      param=ScanBamParam(what=what, which=which)
## Case of no reads on the chromosome   
      start=extractFromBam(file=file_,which=which,what="pos")
      if(length(start)==0 ) {
	brin_F[[i]]=list(numeric(),numeric(),numeric())
	brin_R[[i]]=list(numeric(),numeric(),numeric())
	brin_F_junction[[i]]=list(numeric(),numeric(),numeric())
	brin_R_junction[[i]]=list(numeric(),numeric(),numeric())
      }
      else {
	strand=extractFromBam(file=file_,which=which,what="strand")
	flag=extractFromBam(file=file_,which=which,what="flag")
	mpos=extractFromBam(file=file_,which=which,what="mpos")
	cigar=extractFromBam(file=file_,which=which,what="cigar")
	mrnm=extractFromBam(file=file_,which=which,what="mrnm")
	isize=extractFromBam(file=file_,which=which,what="isize")
	
	##pour inverser les deux brins + et - (* = variable temporaire)
	first_read=integer.base.b(flag)[,5]==1
	strand[strand=="+" & !first_read ]="*"
	strand[strand=="-" & !first_read ]="+"
	strand[strand=="*" & !first_read ]="-"

## CIGAR's interpreter
	cigar_interprete=interpreteCIGAR(cigar)
	longueur_cigar=lapply(cigar_interprete,calcule_longueur_cigar)
	junction_cigar=lapply(cigar_interprete,calcule_junction_cigar)
	
	end=start+sapply(1:length(cigar),function(j) longueur_cigar[[cigar[j]]])

## Case of pairend reads
	##on regarde si R1 et R2 sur meme chromosome
	is_on_same_chr=mrnm==nom_chromo	
	is_on_same_chr[is.na(is_on_same_chr)]=FALSE
	
	##faut que la distance R1 et R2 soit < insert_max
	is_paired=is_on_same_chr & abs(isize) <=insert_max_
	is_paired[first_read & strand=="+" & (isize<0 | isize>insert_max_)]=FALSE
	is_paired[!first_read & strand=="+" & (isize>0 | isize < -insert_max_)]=FALSE
	is_paired[first_read & strand=="-" & (isize>0 | isize < -insert_max_)]=FALSE
	is_paired[!first_read & strand=="-" & (isize<0 | isize>insert_max_)]=FALSE
	is_paired[is.na(is_paired)]=FALSE
	
	##on recupere debut et fin fragment selon la position de R1 et R2
	##R1 sur + et R2 sur -
	##R1 = first_read et R2 = !first_read
	debut_fragment_paired_plus=mpos[!first_read & strand=="+" & is_paired]
	fin_fragment_paired_plus=end[!first_read & strand=="+" & is_paired]
	##R1 sur - et R2 sur +
	debut_fragment_paired_moins=mpos[first_read & strand=="-" & is_paired]
	fin_fragment_paired_moins=end[first_read & strand=="-" & is_paired]
	
## Case of single reads
	##si R1 ou R2 a pas mappé
	debut_fragment_singleton_plus=start[!is_paired & strand=="+"] 
	fin_fragment_singleton_plus=end[!is_paired & strand=="+"]
	debut_fragment_singleton_moins=start[!is_paired & strand=="-"] 
	fin_fragment_singleton_moins=end[!is_paired & strand=="-"]

## Fragments
	debut_frag_plus=c(debut_fragment_paired_plus,debut_fragment_singleton_plus)
	fin_frag_plus=c(fin_fragment_paired_plus,fin_fragment_singleton_plus)
	debut_frag_moins=c(debut_fragment_paired_moins,debut_fragment_singleton_moins)
	fin_frag_moins=c(fin_fragment_paired_moins,fin_fragment_singleton_moins)
	brin_F[[i]]=compresse_coordonnees(debut_frag_plus,fin_frag_plus)
	brin_R[[i]]=compresse_coordonnees(debut_frag_moins,fin_frag_moins)

## Junction read
	debut_junction=numeric()
	fin_junction=numeric()
	brin_junction=numeric()
	i_junction=0
	for(j in 1:length(cigar)) {
	  junctions_=junction_cigar[[cigar[j]]]
	  if(length(junctions_)) {
	    for(k in 1:length(junctions_)) {
	      i_junction=i_junction + 1
	      debut_junction[i_junction]=start[j] + junctions_[[k]][1] - 1
	      fin_junction[i_junction]=start[j] + junctions_[[k]][2] - 1
	      brin_junction[i_junction]=as.character(strand[j])
	    }
	  }
	}
	if(i_junction==0) {
	  brin_F_junction[[i]]=list(numeric(),numeric(),numeric())
	  brin_R_junction[[i]]=list(numeric(),numeric(),numeric())
	}
	else {
	  brin_F_junction[[i]]=compresse_coordonnees(debut_junction[brin_junction=="+"],fin_junction[brin_junction=="+"])
	  brin_R_junction[[i]]=compresse_coordonnees(debut_junction[brin_junction=="-"],fin_junction[brin_junction=="-"])
	}

## Flagstat interpreter
	flag_bits=integer.base.b(flag)#remplie les  données stat pour un flag donné

	## flagstat
	## total
	flagstat[1]=flagstat[1] + sum(flag_bits[,2]==0)
	## duplicates
	flagstat[2]=flagstat[2] + sum((flag_bits[,1]==1)&(flag_bits[,2]==0))
	## mapped
	flagstat[3]=flagstat[3] + sum((flag_bits[,9]==0)&(flag_bits[,2]==0))
	## paired
	flagstat[4]=flagstat[4] + sum((flag_bits[,11]==1)&(flag_bits[,2]==0))
	## read1
	flagstat[5]=flagstat[5] + sum((flag_bits[,5]==1)&(flag_bits[,2]==0))
	## read2
	flagstat[6]=flagstat[6] + sum((flag_bits[,4]==1)&(flag_bits[,2]==0))
	## iself and mate mapped
	flagstat[8]=flagstat[8] + sum((flag_bits[,11]==1)&(flag_bits[,9]==0)&(flag_bits[,8]==0)&(flag_bits[,2]==0))
	## singletons
	flagstat[9]=flagstat[9] + sum((flag_bits[,8]==1)&(flag_bits[,2]==0))
	## QC-failed
	flagstat[11]=flagstat[11] + sum(flag_bits[,2]==1)
	## flagstat
	## mate on a different chr
	flagstat[10]=flagstat[10] + sum((!is_on_same_chr)&(flag_bits[,11]==1)&(flag_bits[,9]==0)&(flag_bits[,8]==0)&(flag_bits[,2]==0))
	## flagstat
	## properly paired
	flagstat[7]=flagstat[7] + sum(is_paired)
      }
    }
  }
  
## Data storing
  names(brin_F)=noms_chromosomes
  names(brin_R)=noms_chromosomes
  names(brin_F_junction)=noms_chromosomes
  names(brin_R_junction)=noms_chromosomes

  bamHandler=list()
  if(libraryType_[1]=="inverse") {
    bamHandler[[1]]=brin_R
    bamHandler[[2]]=brin_F
  }
  else {
    bamHandler[[1]]=brin_F
    bamHandler[[2]]=brin_R
  }
  bamHandler[[3]]=longueur_chromosomes
  bamHandler[[4]]=flagstat
  bamHandler[[5]]=stranded_

  if(libraryType_[1]=="inverse") {
    bamHandler[[6]]=brin_R_junction
    bamHandler[[7]]=brin_F_junction
  }
  else {
    bamHandler[[6]]=brin_F_junction
    bamHandler[[7]]=brin_R_junction
  }
  bamHandler[[8]]=FALSE
  if(!is.null(normalized_)) {
    for( i in pas) {
      for(j in 1:nombre_chromosomes) {
	bamHandler[[i]][[j]][[4]]=normalized_*bamHandler[[i]][[j]][[3]]
      }
    }
    bamHandler[[8]]=TRUE
  }
    
  names(bamHandler)=c("F","R","chrLength","flagstat","stranded","junctions_F","junctions_R","norm")
  if((is.null(chrName_))&(from_==1)&(is.null(to_))) {
    if(is.null(fileNameRData_)|is.na(fileNameRData_)) {
      save(bamHandler,file=RDataFileName(file_))
    } 
    else {
      save(bamHandler,file=fileNameRData_)
    }
  }
  return(bamHandler)
}

chercheBamHandlerDansFichierRData=function(fichierRData) {
  tmp=try(load(fichierRData),TRUE)
  if(class(tmp)=="try-error") {
    return(tmp)
  }
  else {
    return(try(bamHandler,TRUE))
  }
}

## Public Function
readBam=function(file,insert_max=2000,stranded=TRUE,reload=FALSE,ncore=1,libraryType=c("standard","inverse"),normalized=NULL,chrName=NULL,from=1,to=NULL) {
  if(length(file)==1) { ##si un seul fichier
    nom_fichier_RData=RDataFileName(file)
    if(!file.exists(nom_fichier_RData)|reload) {
      return(try(readBam_(file_=file,insert_max_=insert_max,stranded_=stranded,ncore_=ncore,libraryType_=libraryType,normalized_=normalized,chrName_=chrName,from_=from,to_=to),TRUE))
    }
    else {
      load(nom_fichier_RData)
      return(bamHandler)
    }
  }
  else { ##si plusieurs fichiers
    chargeBamsDansFichiersRData=function(oneFile) {
      nom_fichier_RData=RDataFileName(oneFile)
      if(!file.exists(nom_fichier_RData)|reload) {
	bamHandler=try(readBam_(oneFile,insert_max_=insert_max,stranded_=stranded,libraryType_=libraryType,normalized_=normalized,chrName_=chrName,from_=from,to_=to),TRUE)
      }
      return(nom_fichier_RData)
    }
    if(ncore>1) {    
      require("multicore")
      fichiersRData=mclapply(file,chargeBamsDansFichiersRData,mc.cores=ncore,mc.preschedule=FALSE)
    }
    else {
      fichiersRData=lapply(file,chargeBamsDansFichiersRData)
    }
    gc()
    return(lapply(fichiersRData,chercheBamHandlerDansFichierRData))
  }
}

## Returns an adjusted bamHandler object
ajustBam=function(bamHandler,coeff=1) {
  if(class(bamHandler)=="try-error") {
    return(bamHandler)
  }
  else {
    bamHandler_=bamHandler
    L=length(bamHandler_$F)
    if(!bamHandler_$norm) {
      for(i in 1:L) {
	bamHandler_$F[[i]][[4]]=bamHandler_$F[[i]][[3]]*coeff
	bamHandler_$R[[i]][[4]]=bamHandler_$R[[i]][[3]]*coeff
	bamHandler_$junctions_F[[i]][[4]]=bamHandler_$junctions_F[[i]][[3]]*coeff
	bamHandler_$junctions_R[[i]][[4]]=bamHandler_$junctions_R[[i]][[3]]*coeff
      }
      bamHandler_$norm=TRUE
    }
    else {
      for(i in 1:L) {
	bamHandler_$F[[i]][[4]]=bamHandler_$F[[i]][[4]]*coeff
	bamHandler_$R[[i]][[4]]=bamHandler_$R[[i]][[4]]*coeff
	bamHandler_$junctions_F[[i]][[4]]=bamHandler_$junctions_F[[i]][[4]]*coeff
	bamHandler_$junctions_R[[i]][[4]]=bamHandler_$junctions_R[[i]][[4]]*coeff
      }
    }
    return(bamHandler_)
  }
}

## Returns an ajusted bamHandler list
normalizeBams=function(data,poids,which=1:length(poids),relative=TRUE) {
  N=length(data)
  j=0
  retour=list()
  for(i in 1:N) {
    if(i %in% which) {
      j=j + 1
      if(relative) {
	coeff=min(poids)/poids[j]
      }
      else {
	coeff=poids[j]
      }
      retour[[i]]=ajustBam(data[[i]],coeff)
    }
    else {
      retour[[i]]=data[[i]]
    }
  }
  return(retour)
}

## Returns the sum of two bamHandler objects
addBam=function(bamHandler1,bamHandler2) {
  if(class(bamHandler1)=="try-error"|class(bamHandler2)=="try-error") {
    if(class(bamHandler1)=="try-error") {
      return(bamHandler2)
    }
    else {
      return(bamHandler1)
    }
  }
  else {
    brin_F=list()
    junctions_brin_F=list()
    brin_R=list()
    junctions_brin_R=list()
    L=length(bamHandler1$F)
    bamHandler=list()
    if(bamHandler1$norm !=bamHandler2$norm ) {
      warning(expr="Two different bam files(normalized and non normalized)!!",immediate.=TRUE)
    }
    for(i in 1:L) {
      brin_F[[i]]=list(c(bamHandler1$F[[i]][[1]],bamHandler2$F[[i]][[1]]),c(bamHandler1$F[[i]][[2]],bamHandler2$F[[i]][[2]]),c(bamHandler1$F[[i]][[3]],bamHandler2$F[[i]][[3]]))
      junctions_brin_F[[i]]=list(c(bamHandler1$junctions_F[[i]][[1]],bamHandler2$junctions_F[[i]][[1]]),c(bamHandler1$junctions_F[[i]][[2]],bamHandler2$junctions_F[[i]][[2]]),c(bamHandler1$junctions_F[[i]][[3]],bamHandler2$junctions_F[[i]][[3]]))
      brin_R[[i]]=list(c(bamHandler1$R[[i]][[1]],bamHandler2$R[[i]][[1]]),c(bamHandler1$R[[i]][[2]],bamHandler2$R[[i]][[2]]),c(bamHandler1$R[[i]][[3]],bamHandler2$R[[i]][[3]]))
      junctions_brin_R[[i]]=list(c(bamHandler1$junctions_R[[i]][[1]],bamHandler2$junctions_R[[i]][[1]]),c(bamHandler1$junctions_R[[i]][[2]],bamHandler2$junctions_R[[i]][[2]]),c(bamHandler1$junctions_R[[i]][[3]],bamHandler2$junctions_R[[i]][[3]]))
      if(bamHandler1$norm & bamHandler2$norm) {
	brin_F[[i]][[4]]=c(bamHandler1$F[[i]][[4]],bamHandler2$F[[i]][[4]])
	junctions_brin_F[[i]][[4]]=c(bamHandler1$junctions_F[[i]][[4]],bamHandler2$junctions_F[[i]][[4]])
	brin_R[[i]][[4]]=c(bamHandler1$R[[i]][[4]],bamHandler2$R[[i]][[4]])
	junctions_brin_R[[i]][[4]]=c(bamHandler1$junctions_R[[i]][[4]],bamHandler2$junctions_R[[i]][[4]])
      }
    }
    names(brin_F)=names(bamHandler1$F)
    names(brin_R)=names(bamHandler2$R)
    names(junctions_brin_F)=names(bamHandler1$junctions_F)
    names(junctions_brin_R)=names(bamHandler2$junctions_R)
    bamHandler[[1]]=brin_F
    bamHandler[[2]]=brin_R
    bamHandler[[3]]=bamHandler1[[3]]
    bamHandler[[4]]=bamHandler1[[4]] + bamHandler2[[4]]
    bamHandler[[5]]=bamHandler1[[5]] & bamHandler2[[5]]
    bamHandler[[6]]=junctions_brin_F
    bamHandler[[7]]=junctions_brin_R 
    bamHandler[[8]]=bamHandler1$norm & bamHandler2$norm
    names(bamHandler)=c("F","R","chrLength","flagstat","stranded","junctions_F","junctions_R","norm")
    return(bamHandler)
  }
}

## Returns the sum of a list of bamHandler objects
sumBam=function(...) {
  args=list(...)
  retour=args[[1]]
  if(length(args)>1) {
    for(i in 2:length(args)) {
      retour=addBam(retour,args[[i]])
    }
  }
  return(retour)
}

## Returns the mean of bamHandler objects
meanBam=function(...) {
  args=list(...)
  if(length(args) == 1 )
    args <- args[[1]]
  bamsReduced=normalizeBams(args,rep(1,length(args))/length(args),relative=FALSE)
  retour=bamsReduced[[1]]
  if(length(bamsReduced)>1) {
    for(i in 2:length(bamsReduced)) {
      retour=addBam(retour,bamsReduced[[i]])
    }
  }
  return(retour)
}

## Intern function for 
totalReads=function(bamHandler) {
  return(bamHandler[[4]][1]+bamHandler[[4]][11])
}

## Intern function for readGff function
my.read.lines2=function(fname) {
  s=file.info( fname )$size 
  buf=readChar( fname, s, useBytes=T)
  strsplit( buf,"\n",fixed=T,useBytes=T)[[1]]
}

## Extracts the annotation infos from Gff file
readGff=function(file_in, from=1, to=Inf, chr=NULL, infoName=c("ID","Name","Parent","gene","Alias","orf_classification","Ontology_term","Note","GO")) {
  tmp=try(my.read.lines2(file_in))
  if(!is.null(chr)) {
    tmp1=grep(paste(chr,"\t",sep=""), tmp, value=TRUE,useBytes=T)
  }
  else {
    tmp1=tmp
  }
  N=length(tmp1)
  Chr=array()
  Start=array()
  Stop=array()
  Strand=array()
  Type=array()
  info=list()
  for(i in 1:length(infoName)) info[[i]]=array()
  names(info)=infoName
  j=1
  for (i in 1:N) {
    if(substr(tmp1[i],1,1)!="#") {
      line_split=unlist(strsplit(tmp1[i],"\t",fixed=T,useBytes=T))
      if((as.integer(line_split[4])<=to) & (as.integer(line_split[5])>=from)) {
	Chr[j]=line_split[1]
	Start[j]=as.integer(line_split[4])
	Stop[j]=as.integer(line_split[5])
	Strand[j]=line_split[7]
	Type[j]=line_split[3]
	ninth=unlist(strsplit(line_split[9],";",fixed=T,useBytes=T))
	element_ninth_empty=rep(TRUE,length(infoName))
	for(element_ninth in ninth) {
	  element_ninth_split=unlist(strsplit(element_ninth,"=",fixed=T,useBytes=T))
	  if(length(element_ninth_split)==2) {
	    if(element_ninth_split[1] %in% infoName) {
	      info[[element_ninth_split[1]]][j]=element_ninth_split[2]
	      element_ninth_empty[infoName==element_ninth_split[1]]=FALSE
	    }
	  }
	}
	for(infoName_ in infoName[element_ninth_empty]) {
	  info[[infoName_]][j]="."
	}
	j=j+1
      }
    }
  } 
  retour=data.frame(Chr,Type,Start,Stop,Strand,info,stringsAsFactors=FALSE)
  return(retour)
}

## Counts the Reads for each transcript for one bamHandler object
countReads=function(bamHandler,annot,label="",normalized_=FALSE) {
  N=dim(annot)[1]
  densite=rep(0,N)
  nombrereads=rep(0,N)
  if(normalized_) {
    if(bamHandler$norm) {
      j=4
    }else {
      j=3
    }
  }else {
    j=3
  }
  if(bamHandler$stranded) {
    for(i in 1:N) {
      chr_=annot$Chr[i]
      if(annot$Strand[i]=="+") {
	if(!is.null(bamHandler$F[[chr_]])) {
	  tmp=sum(bamHandler$F[[chr_]][[j]][(bamHandler$F[[chr_]][[1]]<=annot$Stop[i])&(bamHandler$F[[chr_]][[2]]>=annot$Start[i])])
	}
	else {
	  tmp=0
	}
	densite[i]=tmp/(annot$Stop[i]-annot$Start[i]+1)
	nombrereads[i]=tmp
      }
      if(annot$Strand[i]=="-") {
	if(!is.null(bamHandler$R[[chr_]])) {   
	  tmp=sum(bamHandler$R[[chr_]][[j]][(bamHandler$R[[chr_]][[1]]<=annot$Stop[i])&(bamHandler$R[[chr_]][[2]]>=annot$Start[i])])
	}
	else {
	  tmp=0
	}
	densite[i]=tmp/(annot$Stop[i]-annot$Start[i]+1)
	nombrereads[i]=tmp
      }
      if(annot$Strand[i]==".") {
	if(!is.null(bamHandler$F[[chr_]])) {
	  tmpF=sum(bamHandler$F[[chr_]][[j]][(bamHandler$F[[chr_]][[1]]<=annot$Stop[i])&(bamHandler$F[[chr_]][[2]]>=annot$Start[i])])
	}
	else {
	  tmpF=0
	}
	if(!is.null(bamHandler$R[[chr_]])) {
	  tmpR=sum(bamHandler$R[[chr_]][[j]][(bamHandler$R[[chr_]][[1]]<=annot$Stop[i])&(bamHandler$R[[chr_]][[2]]>=annot$Start[i])])
	}
	else {
	  tmpR=0
	}
	densite[i]=(tmpF+tmpR)/(annot$Stop[i]-annot$Start[i]+1)
	nombrereads[i]=tmpF+tmpR
      }
    }
  }
  else {
    for(i in 1:N) {
      chr_=annot$Chr[i]
      if(!is.null(bamHandler$F[[chr_]])) {
	tmpF=sum(bamHandler$F[[chr_]][[j]][(bamHandler$F[[chr_]][[1]]<=annot$Stop[i])&(bamHandler$F[[chr_]][[2]]>=annot$Start[i])])
      }
      else {
	tmpF=0
      }
      if(!is.null(bamHandler$R[[chr_]])) {
	tmpR=sum(bamHandler$R[[chr_]][[j]][(bamHandler$R[[chr_]][[1]]<=annot$Stop[i])&(bamHandler$R[[chr_]][[2]]>=annot$Start[i])])
      }
      else {
	tmpR=0
      }
      densite[i]=(tmpF+tmpR)/(annot$Stop[i]-annot$Start[i]+1)
      nombrereads[i]=tmpF+tmpR
    }
  }     
  retour=data.frame(densite,nombrereads)
  if(label=="") {
    names(retour)=c("densities","readcount")
  }
  else {
    names(retour)=paste(label,c("densities","readcount"),sep=" ")
  }

  return(retour)
}

## Counts the Reads for each transcript for different bamHandler objects
countReadsBams=function(data,description_data,annot,ncore=1,which=1:length(data),normalized=FALSE) {
  doit=function(i) {
    return(countReads(data[[i]],annot,label=description_data[i],normalized_=normalized))
  }
  if(ncore==1) {
    list_retour=lapply(which,doit)
  }
  else {
    list_retour=mclapply(which,doit,mc.cores=ncore,mc.preschedule=FALSE)
  }
  return(listOfDataFrame2DataFrame(list_retour,vertical=FALSE))
}

##### Other functions

# function to read output files from bedtools intersect (intersect of gff files)
# input : bedtools intersect output 
# output : data.frame, with only ID filed kept from last gff column 

readIntersect <- function(file){
  inter <- read.csv(file, header = F, sep = "\t", stringsAsFactors = F)
  inter$V9 <- gsub("ID=", "", unlist(lapply(inter$V9, function(l){
    tmp <- strsplit(as.character(l), ";")
    return(tmp[[1]][1])})))
  inter$V18 <- gsub("ID=", "", unlist(lapply(inter$V18, function(l){
    tmp <- strsplit(as.character(l), ";")
    return(tmp[[1]][1])})))
  return(inter)
}

# function to write in gff file

writeGff <- function(df, out){
  write.table(cbind.data.frame(chr=df$Chr,
                               source = rep(".", nrow(df)),
                               Type = df$Type,
                               Start = df$Start,
                               Stop = df$Stop,
                               rep(".", nrow(df)),
                               Strand = df$Strand,
                               rep(".", nrow(df)),
                               ID = paste0("ID=", df$ID)), 
              file = out, quote = F, sep = "\t", row.names = F, col.names = F)
}
