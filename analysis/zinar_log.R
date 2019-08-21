#  *************************************************************************************************************
#  *  Initially written by Marc Descrimes                                                                      *
#  *                                                                                                           *
#  *  script to segment signal :                                                                               *
#  *  Coverage is computed at the nucleotide level input                                                       *
#  *  for each window size and threshold given, a gff file is produced with coordinates of segments for which  *
#  *  the mean of signal (log2) over the window is above threshold                                             *             
#  *************************************************************************************************************

options(scipen = 7)
## Get Arguments
set = function(nomVariable,valeur) {
  if(exists(nomVariable)) {
    L=length(get(nomVariable))
  }
  else {
    L=0
    commande = paste(nomVariable,"<<-numeric()",sep="")
    expr=try(parse(text=commande),TRUE)
    eval(expr)
  }
  if(is.character(valeur)) {
    commande = paste(nomVariable,"[",L+1,"]<<-\"",valeur,"\"",sep="")
  }
  else {
    commande = paste(nomVariable,"[",L+1,"]<<-",valeur,sep="")
  }
  expr=try(parse(text=commande),TRUE)
  eval(expr)
}

nomValParam = c()
defautParam = list()
attribParam = c()
descriParam = c()
typeParam = c()
n_param = 0

newParam = function(nom,def,att,desc,type) {
  n_param <<- n_param+1
  nomValParam[n_param] <<- nom
  defautParam[n_param] <<- def
  attribParam[n_param] <<- att
  descriParam[n_param] <<- desc
  typeParam[n_param] <<- type
}

#paramètres
newParam("PARAM_L",60,"-L","Taille de la fenêtre glissante","integer")
newParam("PARAM_s",0.18,"-s","Seuil de détection des transcrits","numeric")
newParam("PARAM_nCore",1,"-c","Nombre de coeur","integer")
newParam("PARAM_prefix","zinar","-n","Préfixe pour le noms des transcrits","character")
newParam("PARAM_src","zinar","--source","Nom de la source des données","character")
newParam("PARAM_imin",50,"--imin","Taille d'insert minimum des paired read","integer")
newParam("PARAM_imax",3000,"--imax","Taille d'insert maximum des paired read","integer")
newParam("PARAM_limIntron",2000,"-I","Taille d'intron maximum","integer")
newParam("PARAM_repertoire","zinar","-o","Chemin du repertoire de sortie","character")
newParam("PARAM_vary",0,"-v","Système vary©","integer")
newParam("PARAM_taille_region",500000,"-t","Taille max d\'une region","integer")
newParam("PARAM_insert_max",2000,"-i","Taille max d\'un insert","integer")
newParam("PARAM_library_type","normal","-l","Type de library","character")

printHelp = function() {
  cat("Zinar v0.01
Commande : zinar.R [options] targetFile
Options :
")
  for(i in 1:n_param) {
    cat(paste(attribParam[i]," :\t",descriParam[i]," (Defaut : ",sep=""))
    cat(paste(defautParam[[i]],")\n",sep=""))
  }
  q("no")
}

#Attributions des valeurs utilisateurs
args = commandArgs(trailingOnly=TRUE)

{
nbr_arguments = length(args)
if(nbr_arguments == 0) {
  printHelp()
}
else {
  bamFile = args[nbr_arguments]
  if(nbr_arguments>1) {
    mode=0
    for(i in 1:(nbr_arguments-1)) {
      if(mode==0) {
	quelParam = which(args[i]==attribParam)
	if(length(quelParam)==0) {
	  cat(paste(args[i]," n'est pas un argument valide\n",sep=""))
	  printHelp()
	}
	else {
	  mode=1
	  paramAttrib = nomValParam[quelParam]
	}
      }
      else {
	set(paramAttrib,args[i])
	mode=0
      }
    }
    if(mode==1) {
      cat("Mauvais appel du programme\n")
      printHelp()
    }
  }
}
}

#Attributions des valeurs par défaut
for(i in 1:n_param) {
  if(!exists(nomValParam[i])) {
    set(nomValParam[i],defautParam[[i]])
  }
}

#options(warn=-1)
#Retypage des arguments
convert = function(x,type) {
  switch( type,
	  numeric = as.numeric(x),
	  integer = as.integer(x),
	  character = as.character(x),
	  boolean = x
	)
}
for(i in 1:n_param) {
  bonneVal = convert(get(nomValParam[i]),typeParam[i])
  if(sum(is.na(bonneVal))!=0) {
    cat(paste("Erreur : ",nomValParam[i]," doit être de type ",typeParam[i],"\n",sep=""))
  }
  else {
    commande = paste(nomValParam[i],"<-bonneVal",sep="")
    expr=try(parse(text=commande),TRUE)
    eval(expr)
  }
}

if(PARAM_vary==1) {
  if(length(PARAM_L)==3) {
    PARAM_L = seq(from=PARAM_L[1],to=PARAM_L[2],by=PARAM_L[3])
  }
  if(length(PARAM_s)==3) {
    PARAM_s = seq(from=PARAM_s[1],to=PARAM_s[2],by=PARAM_s[3])
  }  
}

#debug
for(i in 1:n_param) {
  cat(paste(nomValParam[i],"\t",get(nomValParam[i]),"\t",class(get(nomValParam[i])),"\n",sep=""))
}
cat(paste("bamFile\t",bamFile,"\n",sep=""))

##MAIN
#Fonction charge et quitte si le package n'est pas installé
charge = function(truc) {
  x=suppressPackageStartupMessages(require(truc,quietly=TRUE,character.only=TRUE,warn.conflicts=FALSE))
  if(x==FALSE) {
    warning(paste("Il faut installer le package ",truc,"\n",sep=""))
    q("no")
  }
}
#Traitement simply
charge("Rsamtools")

#fonction conversion binaire
integer.base.b = function(x, b=2){
  xi <- as.integer(x)
  if(any(is.na(xi) | ((x-xi)!=0)))
    print(list(ERROR="x not integer", x=x))
  N <- length(x)
  xMax <- max(x)
  ndigits <- 32
  Base.b <- array(NA, dim=c(N, ndigits))
  for(i in 1:ndigits){#i <- 1
    Base.b[, ndigits-i+1] <- (x %% b)
    x <- (x %/% b)
  }
  Base.b
}

#fonction qui convertit une liste de data frame vers un data frame
listOfDataFrame2DataFrame = function(x) {
  if(length(x)==1) {
    x[[1]]
  }
  else {
    ncol = dim(x[[1]])[2]
    vect=unlist(lapply(x,t))
    retour = as.data.frame(t(array(vect,dim=c(ncol,length(vect)/ncol))),stringsAsFactors=FALSE)
    names(retour)=names(x[[1]])
    retour
  }
}

#Récupération info bam
genome_info=scanBamHeader(bamFile)[[1]]$targets
noms_chromosomes = names(genome_info)
longueur_chromosomes = as.integer(genome_info)
nombre_chromosomes = length(noms_chromosomes)

n_region_par_chromosome=round(longueur_chromosomes/PARAM_taille_region)
n_region_par_chromosome[n_region_par_chromosome==0]=1

# print(cbind(longueur_chromosomes,n_region_par_chromosome))

chromo_region=c()
start_region=c()
stop_region=c()

for(i in 1:nombre_chromosomes) {
  if(n_region_par_chromosome[i]==1) {
    chromo_region=c(chromo_region,noms_chromosomes[i])
    start_region=c(start_region,1)
    stop_region=c(stop_region,longueur_chromosomes[i])
  }
  else {
    coupures=ceiling((1:n_region_par_chromosome[i])*(longueur_chromosomes[i]/(n_region_par_chromosome[i])))
    chromo_region=c(chromo_region,noms_chromosomes[i])
    start_region=c(start_region,1)
    stop_region=c(stop_region,coupures[1])
    for(j in 2:length(coupures)) {
      chromo_region=c(chromo_region,noms_chromosomes[i])
      start_region=c(start_region,coupures[j-1]+1)
      stop_region=c(stop_region,coupures[j])
    }
  }
}

nombre_regions=length(chromo_region)
# nombre_regions=120

# print(data.frame(chromo_region,start_region,stop_region))

check_bam = function(i) {
  retour=list()

  nom_chromo = chromo_region[i]
  start_region_ = start_region[i]
  stop_region_ = stop_region[i]

  commande = paste("RangesList(`",nom_chromo,"` = IRanges(",start_region_,",",stop_region_,"))",sep="")
  
  expr=try(parse(text=commande),TRUE)

  if(class(expr)=="expression") {
    which = eval(expr)

    flag = scanBam(bamFile, param=ScanBamParam(which=which,what="flag"))[[1]][[1]]
    if(length(flag)==0){
      retour[[1]] = 0
      retour[[2]] = 0
      retour[[3]] = NA
      retour[[4]] = NA
      retour[[5]] = 0      
    }
    else {
      qwidth = scanBam(bamFile, param=ScanBamParam(which=which,what="qwidth"))[[1]][[1]]

      flag_bits = integer.base.b(flag)

      n_reads = dim(flag_bits)[1]
      n_bits = dim(flag_bits)[2]

      multiple = flag_bits[,n_bits]
      reverse = flag_bits[,n_bits-4]
      first = flag_bits[,n_bits-6]

      doit = function(j) {
	if(multiple[j] == 0) {
	  if(reverse[j] == 0) {
	    0
	  }
	  else {
	    1
	  }
	}
	else {
	  if(reverse[j] == 0) {
	    if(first[j] == 0) {
	      1
	    }
	    else {
	      0
	    }
	  }
	  else {
	    if(first[j] == 0) {
	      0
	    }
	    else {
	      1
	    }
	  }
	}
      }
      brin = unlist(lapply(1:n_reads,doit))
      
      retour[[1]] = n_reads
      retour[[2]] = sum(multiple)
      retour[[3]] = unique(qwidth[first==0])
      retour[[4]] = unique(qwidth[first==1])
      retour[[5]] = sum(qwidth)
    }
  }
  else {
    print("nope")
  }

  return(retour)
}

cat("\nContrôle des reads\n")
{
  if(PARAM_nCore[1]>1) {
    charge("multicore")
    controle=mclapply(1:nombre_regions,check_bam,mc.cores=PARAM_nCore[1],mc.preschedule=FALSE)
  }
  else {
    controle=lapply(1:nombre_regions,check_bam)
  }
}

toto=unlist(lapply(controle,function(x) {if(class(x)=="try-error") {print(x);return(1)} else {return(0)}}))

nombre_reads_total = 0
nombre_reads_paired = 0
nombre_nucleotides_lus = 0

for(i in 1:nombre_regions) {
  if(length(controle[[i]])!=0) {
    nombre_reads_total = nombre_reads_total + controle[[i]][[1]]
    nombre_reads_paired = nombre_reads_paired + controle[[i]][[2]]
    nombre_nucleotides_lus = nombre_nucleotides_lus + controle[[i]][[5]]
  }
}

longueur_genome = sum(as.numeric(longueur_chromosomes))
couverture = nombre_nucleotides_lus/longueur_genome

{
  if(nombre_reads_paired == 0) {
    type_library = "single"
    longueurs_reads = c()
    for(i in 1:nombre_chromosomes) {
      if(length(controle[[i]])!=0) {
	longueurs_reads = unique(c(longueurs_reads,controle[[i]][[3]]))
      }
    }
    longueurs_reads_homogene = length(longueurs_reads)==1
    cat("Données Single-Read\n")
    cat(paste("Reads length : ",paste(longueurs_reads,collapse=" "),"\n"))
  }
  else {
    type_library = "paired"
    longueurs_reads = c()
    longueurs_reads_mate = c()
    for(i in 1:nombre_chromosomes) {
      if(length(controle[[i]])!=0) {
	longueurs_reads = unique(c(longueurs_reads,controle[[i]][[4]]))
	longueurs_reads_mate = unique(c(longueurs_reads_mate,controle[[i]][[3]]))
      }
    }
    cat("Données Paired-End\n")
    cat(paste("Reads length first mate  : ",paste(longueurs_reads,collapse=" "),"\n"))
    cat(paste("Reads length second mate : ",paste(longueurs_reads_mate,collapse=" "),"\n"))
  }
  # if(couverture>=10) {
  #   cat(paste("Couverture : ",round(couverture),"X\n",sep=""))
  # }
  # else {
  #   cat(paste("Couverture : ",round(couverture,2),"X\n",sep=""))  
  # }
}

# save.image("zinar_test_1.RData")
# stop()

#Fonctions utiles
extractFromBam = function(bamFile,which,what) {
  return(scanBam(bamFile, param=ScanBamParam(which=which,what=what))[[1]][[1]])
}

interpreteCIGAR = function(cigar) {
  cigar_un = strsplit(unique(cigar),split="")
  n_cigar_un = length(cigar_un)
  taille_cigar = list()
  analise_cigar = function(cigar_) {
    cigar_sortie = list()
    acc = ""
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
  cigar_interprete = lapply(cigar_un,analise_cigar)
  names(cigar_interprete) = unique(cigar)

  return(cigar_interprete)
}

calcule_longueur_cigar = function(cigar_) {
  lon = 0
  N=length(cigar_)
  for(j in seq(2,N,2)) {
    if(cigar_[[j]]!="I") {
      lon=lon+cigar_[[j-1]]
    }
  }
  return(lon)
}

segment = function(signal, chrname, chrstrand,start_region_) {
# x : trucs à segmenter
# PARAM_L : zone de lissage (60)
# PARAM_s : seuil de segment potentiel (0.09)
# PARAM_d : distance de fusion potentielle (250)
# PARAM_pval : p-value du test de Fisher (0.01)

  n_segmentation = 0
  reponse = list()

  x = signal
  x_log = log2(1+signal)
  clean = function(X) {
    return(X[!is.na(X)])
  }
  long_x = length(x)

  for(i_L in 1:length(PARAM_L)) {
    x_lisse=numeric()
    if(length(x)>0) {
      x_lisse = numeric(long_x)
      if(PARAM_L[i_L]>0) {
	for(i_nucl in 1:long_x) {
	  x_lisse[i_nucl]=mean(x_log[max(1,i_nucl-PARAM_L[i_L]):min(long_x,i_nucl+PARAM_L[i_L])])
	}
      }
      else {
	x_lisse=x_log
      }
    }
    for(i_s in 1:length(PARAM_s)) {
      if(length(x_lisse)>0) {

	x_seuil = c(FALSE,(x_lisse>=log2(1+PARAM_s[i_s])),FALSE)
	x_point = x_seuil[2:(long_x+2)]-x_seuil[1:(long_x+1)]
	x_start = (1:(long_x+1))[x_point==1]
	x_stop = (1:(long_x+1))[x_point==-1]-1

	n_seg = length(x_start)
	if(n_seg>0) {
	  for(i_nucl in 1:n_seg) {
	    seg_start = x_start[i_nucl]
	    seg_stop = x_stop[i_nucl]
	    tmp = (seg_start:seg_stop)[x[seg_start:seg_stop]>0]
	    if(length(tmp)>1) {
	      x_start[i_nucl] = tmp[1]
	      x_stop[i_nucl] = tmp[length(tmp)]
	    }
	    else{
	      x_start[i_nucl] = NA
	      x_stop[i_nucl] = NA
	    }
	  }
	  x_start = clean(x_start)
	  x_stop = clean(x_stop)
	}

	n_seg = length(x_start)
	if(n_seg>0) {
	  chr = rep(chrname,n_seg)
	  strand = rep(chrstrand,n_seg)
	  
	  reponse_ = data.frame(x_start+start_region_-1,x_stop+start_region_-1,strand,chr, stringsAsFactors = FALSE)
	}
	else {
	  reponse_ = data.frame(NA,NA,chrstrand,chrname, stringsAsFactors = FALSE)
	}
      }
      else {
	reponse_ = data.frame(NA,NA,chrstrand,chrname, stringsAsFactors = FALSE)
      }
      names(reponse_) = c("Start","Stop","Strand","Chr")
      n_segmentation = n_segmentation + 1
      reponse[[n_segmentation]] = reponse_
    }
  }

  return(reponse)
}

assemble_frag = function(brin_frag,debut_frag,fin_frag,nom_chromo,lon_chromo,start_region_,stop_region_) {
  
  taille_region_=stop_region_-start_region_+1
  
  signal_F = numeric(taille_region_)
  signal_R = numeric(taille_region_)
  
  for(i_frag in 1:length(brin_frag)) {
    A = max(1,debut_frag[i_frag]-start_region_+1)
    B = min(fin_frag[i_frag]-start_region_+1,taille_region_)
    if(brin_frag[i_frag]==0) {
      signal_F[A:B] = signal_F[A:B]+1
    }
    else {
      signal_R[A:B] = signal_R[A:B]+1
    }
  }
  
  segmentation_F = segment(signal_F, nom_chromo, "+",start_region_)
  segmentation_R = segment(signal_R, nom_chromo, "-",start_region_)
  return(list(segmentation_F,segmentation_R))
}

#Segmentation simpliste single read
seg_single_read = function(i) {
  retour=NA

  nom_chromo = chromo_region[i]
  start_region_ = start_region[i]
  stop_region_ = stop_region[i]

  commande = paste("RangesList(`",nom_chromo,"` = IRanges(",start_region_,",",stop_region_,"))",sep="")
  
  expr=try(parse(text=commande),TRUE)

  if(class(expr)=="expression") {
    which = eval(expr)

    pos = extractFromBam(bamFile,which,"pos")
    if(length(pos)==0) {
      retour=list()
      for(i_seg in 1:n_segmentation) {
	retour[[i_seg]]=list(data.frame(NA,NA,"+",nom_chromo, stringsAsFactors = FALSE),data.frame(NA,NA,"-",nom_chromo, stringsAsFactors = FALSE))
      }
    }
    else {
      flag = extractFromBam(bamFile,which,"flag")
      cigar = extractFromBam(bamFile,which,"cigar")

      flag_bits = integer.base.b(flag)

      n_reads = dim(flag_bits)[1]
      n_bits = dim(flag_bits)[2]

      brin = flag_bits[,n_bits-4]

      #interpréteur de CIGAR
      cigar_interprete = interpreteCIGAR(cigar)
      longueur_cigar = lapply(cigar_interprete,calcule_longueur_cigar)

      #Déterminer les coordonnées de fin de fragments
      debut_frag=pos
      fin_frag=numeric()
      doit = function(j) {return(pos[j] + longueur_cigar[[cigar[j]]])}
      fin_frag = unlist(lapply(1:n_reads,doit))

      #SEGMENTATIOOOOOOOOOON !!!
      retour=assemble_frag(brin,debut_frag,fin_frag,nom_chromo,lon_chromo,start_region_,stop_region_)
    }
  }

  return(retour)
}

#Segmentation simpliste paired end
seg_paired_end = function(i) {
  print(i)
  retour=NA

  nom_chromo = chromo_region[i]
  start_region_ = start_region[i]
  stop_region_ = stop_region[i]

  commande = paste("RangesList(`",nom_chromo,"` = IRanges(",start_region_,",",stop_region_,"))",sep="")
  
  expr=try(parse(text=commande),TRUE)

  if(class(expr)=="expression") {
    which = eval(expr)

    start=extractFromBam(bamFile,which,"pos")
    strand=extractFromBam(bamFile,which,"strand")
    flag=extractFromBam(bamFile,which,"flag")
    mpos=extractFromBam(bamFile,which,"mpos")
    cigar=extractFromBam(bamFile,which,"cigar")
    mrnm=extractFromBam(bamFile,which,"mrnm")
    isize=extractFromBam(bamFile,which,"isize")
    
    flag_bits = integer.base.b(flag)
    first_read=flag_bits[,dim(flag_bits)[2]-6]==1
    strand[strand =="+" & !first_read ]="*"
    strand[strand =="-" & !first_read ] ="+"
    strand[strand =="*" & !first_read ] ="-"

##CIGAR's interpreter
    cigar_interprete = interpreteCIGAR(cigar)
    longueur_cigar = lapply(cigar_interprete,calcule_longueur_cigar)
    
    end=start+sapply(1:length(cigar),function(j) longueur_cigar[[cigar[j]]])

##Case of pairend reads
    is_on_same_chr = mrnm==nom_chromo
    is_on_same_chr[is.na(is_on_same_chr)] = FALSE
    is_paired = is_on_same_chr & abs(isize) <= PARAM_insert_max
    is_paired[first_read & strand=="+" & (isize<0 | isize>PARAM_insert_max)] = FALSE
    is_paired[!first_read & strand=="+" & (isize>0 | isize < -PARAM_insert_max)] = FALSE
    is_paired[first_read & strand=="-" & (isize>0 | isize < -PARAM_insert_max)] = FALSE
    is_paired[!first_read & strand=="-" & (isize<0 | isize>PARAM_insert_max)] = FALSE
    is_paired[is.na(is_paired)] = FALSE
    
    debut_fragment_paired_plus<-mpos[!first_read & strand =="+" & is_paired]
    fin_fragment_paired_plus<-end[!first_read & strand=="+" & is_paired]
    debut_fragment_paired_moins<-mpos[first_read & strand=="-" & is_paired]
    fin_fragment_paired_moins<-end[first_read & strand =="-" & is_paired]
    
##Case of single reads
    debut_fragment_singleton_plus = start[!is_paired & strand =="+"] 
    fin_fragment_singleton_plus = end[!is_paired & strand =="+"]
    debut_fragment_singleton_moins=start[!is_paired & strand =="-"] 
    fin_fragment_singleton_moins= end[!is_paired & strand =="-"]

##Fragments
    debut_frag_plus = c(debut_fragment_paired_plus,debut_fragment_singleton_plus)
    fin_frag_plus  = c(fin_fragment_paired_plus,fin_fragment_singleton_plus)
    debut_frag_moins=c(debut_fragment_paired_moins,debut_fragment_singleton_moins)
    fin_frag_moins = c(fin_fragment_paired_moins,fin_fragment_singleton_moins)

    debut_frag = c(debut_frag_plus,debut_frag_moins)
    fin_frag = c(fin_frag_plus,fin_frag_moins)
    if(PARAM_library_type=="normal") {
      brin_frag = c(rep(0,length(debut_frag_plus)),rep(1,length(debut_frag_moins)))
    }
    else {
      if(PARAM_library_type=="inverse") {
	brin_frag = c(rep(1,length(debut_frag_plus)),rep(0,length(debut_frag_moins)))    
      }
      else {
	brin_frag = rep(1,length(debut_frag_plus)+length(debut_frag_moins))
      }
    }

    orderByStart = order(debut_frag)
    brin_frag=brin_frag[orderByStart]
    debut_frag=debut_frag[orderByStart]
    fin_frag=fin_frag[orderByStart]
    #SEGMENTATIOOOOOOOOOON !!!
    retour=assemble_frag(brin_frag,debut_frag,fin_frag,nom_chromo,lon_chromo,start_region_,stop_region_)
  }

  return(retour)
}

{
  if(type_library == "single") {
    cat("\nSegmentation Single Read\n")
    if(PARAM_nCore[1]>1) {
      charge("multicore")
      seg=mclapply(1:nombre_regions,seg_single_read,mc.cores=PARAM_nCore[1],mc.preschedule=FALSE)
    }
    else {
      seg=lapply(1:nombre_regions,seg_single_read)
    }
  }
  else {
    cat("\nSegmentation Paired End\n")
    if(PARAM_nCore[1]>1) {
      charge("multicore")
      seg=mclapply(1:nombre_regions,seg_paired_end,mc.cores=PARAM_nCore[1],mc.preschedule=FALSE)
    }
    else {
      seg=lapply(1:nombre_regions,seg_paired_end)
    }    
  }
}

save.image("zinar_test.RData")

#Ecrire seg au format gff
system(paste("mkdir -p",PARAM_repertoire[1]))
system(paste("rm -f ",PARAM_repertoire[1],"/info_params.txt",sep=""))
cat("n_segmentation\tPARAM_L\tPARAM_s\n",file=paste(PARAM_repertoire[1],"/info_params.txt",sep=""))
n_segmentation=0
for(i_L in 1:length(PARAM_L)) {
  for(i_s in 1:length(PARAM_s)) {
    n_segmentation=n_segmentation+1
    cat(paste(n_segmentation,"\t",PARAM_L[i_L],"\t",PARAM_s[i_s],"\n",sep=""),file=paste(PARAM_repertoire[1],"/info_params.txt",sep=""),append=TRUE)
  }
}

tmp=data.frame(chromo_region,start_region,stop_region,stringsAsFactors=FALSE)
write.table(tmp,paste(PARAM_repertoire[1],"/coupures_regions.tab",sep=""),quote=FALSE,sep="\t",row.names=FALSE)

cat(paste("Analyses disponibles dans le repertoire ",PARAM_repertoire[1],"\n",sep=""))
for(i_seg in 1:n_segmentation) {
  seg_=NULL
  append=FALSE
  n_label=0
  for(s in 1:2) {
    for(i in 1:length(seg)) {
      if(class(seg[[i]][[1]][[i_seg]])=="data.frame") {
	seg_=rbind(seg_,seg[[i]][[s]][[i_seg]])
      }
      else {
	print(i)
      }
    }
  }
  seg_=seg_[rowSums(is.na(seg_))==0,]
  N=dim(seg_)[1]
  for(j in 2:N) {
    if(seg_$Chr[j-1]==seg_$Chr[j]) {
      if(seg_$Strand[j-1]==seg_$Strand[j]) {
	if(seg_$Stop[j-1]+1==seg_$Start[j]) {
	  seg_$Start[j]=seg_$Start[j-1]
	  seg_$Chr[j-1]=NA
	  seg_$Start[j-1]=NA
	  seg_$Stop[j-1]=NA
	  seg_$Strand[j-1]=NA
	}
      }
    }
  }
  seg_=seg_[!is.na(seg_$Chr),]
  N=dim(seg_)[1]
  chr=seg_$Chr
  start=seg_$Start
  stop=seg_$Stop
  brin=seg_$Strand
  data = rep(PARAM_src[1],N)
  labels = paste("ID=",PARAM_prefix[1],"_",(1:N)+n_label,sep="")
  n_label=n_label+N
  pointless_column=rep(".",N)
  tmp=data.frame(chr,pointless_column,data,start,stop,pointless_column,brin,pointless_column,labels,stringsAsFactors=FALSE)
  write.table(tmp,paste(PARAM_repertoire[1],"/transcrits_",i_seg,".gff",sep=""),append=append,quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
  append=TRUE
  cat(paste(n_label," transcrits trouvés pour la segmentation ",i_seg,".\n",sep=""))
}
