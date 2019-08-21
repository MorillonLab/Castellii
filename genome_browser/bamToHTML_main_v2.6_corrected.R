## titre_page	:	Titre de la page index
## repertoire	:	Nom du répertoire de sortie
## common_files_dir : répertoire des fichiers style, design.css, etc...
## taille_fenetre	: Taille de la fenêtre de visualisation
## decal		Écart entre chaque fenêtre
## style		Data frame : Type(type de transcrit) col(couleur) label(présence des ID sur le graphe)
## annot		Data frame : annotations, sélection ou concaténation de gff (obtenus par readBam)
## nonSGD		Vecteur : type de transcrit à ne pas lier à SGD
## data			list : les différents mapping
## description_data	vecteur : nom des différents mapping
## texte_description	Texte apparaissant sur la page index
## contenu_vue		list : les différentes vues, chaque vue est un vecteur contenant les indices des éléments de data
## description_vue	vecteur : nom des différentes vues
## groupe_maximum	list : indique dans chaque vue les groupe de mapping devant être visualisés avec la même échelle

library("parallel")

cat("Check and initialize variables\n")

{
  if(!exists("data")) stop("data manquant")
  
  if(!exists("noms_chromosome")){
    noms_chromosome=names(data[[1]]$F)
    n_chromosome = length(noms_chromosome)
    tailles_chromosome = data[[1]]$chrLength
  } else {
    n_chromosome = length(noms_chromosome)
    if(!exists("tailles_chromosome")){
      tailles_chromosome=numeric(n_chromosome)
      tailles_all_chromosome <- data[[1]]$chrLength
      noms_all_chromosome <- names(data[[1]]$F)
      for(i in n_chromosome){
        tailles_chromosome[i] <- tailles_all_chromosome[which(noms_all_chromosome == noms_chromosome[i])]
      }
    }
  }
  
  if(!exists("annot")) stop("annot manquant")
  
  if(!exists("index_data")) stop("index_data manquant")
  
  if(!exists("titre_page")) titre_page="sequencing visualisation"
  
  if(!exists("repertoire")) repertoire="out"
  
  if(!exists("taille_fenetre")) taille_fenetre=15000
  
  if(!exists("decal")) decal=5000
  
  if(!exists("style")) {
    Type=unique(annot$Type)
    col=rep(rainbow(length(Type)),2)
    label=1
    shape="box"
    Strand=c(rep("+",length(Type)),rep("-",length(Type)))
    style=data.frame(Type,col,label,shape,Strand,stringsAsFactors = FALSE)
  }
  else {
    style=style[style$Type %in% unique(annot$Type),]
    Type=unique(annot$Type)
    Type_sans_style=Type[!(Type %in% style$Type)]
    if(length(Type_sans_style)>0) {
      col="#FFFFFF"
      label=1
      shape="box"
      Strand=c(rep("+",length(Type_sans_style)),rep("-",length(Type_sans_style)))
      style_complement=data.frame(Type=Type_sans_style,col,label,shape,Strand,stringsAsFactors = FALSE)
      style=rbind(style,style_complement)
    }
  }
  
  if(!exists("nonSGD")) nonSGD = c()
  
  if(!exists("noncount")) noncount = ""
  
  if(!exists("nonform")) nonform = ""
  
  if(!exists("display_name")) display_name = NULL
  
  if(!exists("ncore")) ncore=1
  
  if(!exists("description_data")) description_data=index_data
  
  if(!exists("texte_description")) texte_description=""
  
  if(!exists("contenu_vue")) {
    contenu_vue = list()
    contenu_vue[[1]] = 1:length(index_data)
  }
  
  if(!exists("contenu_vue_analyse_generiq")) {
    contenu_vue_analyse_generiq = contenu_vue
  }
  
  if(!exists("isVisuLog")) {
    isVisuLog=TRUE
  }
  
  if(!exists("isVisuStranded")) {
    isVisuStranded=TRUE
  }
  
  if(!exists("typeVisu")) typeVisu = rep("classic",length(contenu_vue))
  
  if(!exists("classic_plus_color")) {
    classic_plus_color="navyblue"
  }
  
  if(!exists("classic_minus_color")) {
    classic_minus_color="deeppink3"
  }
  
  if(!exists("heatmap_max_color")) {
    heatmap_max_color="#000055"
  }
  
  if(!exists("heatmap_min_color")) {
    heatmap_min_color="#FFFFAA"
  }
  
  if(!exists("heatmap_palette_method")) {
    heatmap_palette_method="hsv"
  }
  
  if(!exists("lines_samples_colors")) {
    lines_samples_colors=c(1,3,4,2)[((0:length(description_data))%%4)+1]
  }
  
  if(!exists("lines_samples_type_line")) {
    lines_samples_type_line=((0:length(description_data))%/%4)+1
  }
  
  if(!exists("smoothLength")) {
    smoothLength=trunc(taille_fenetre/1200)
  }
  
  if(!exists("description_vue")) {
    description_vue=c()
    description_vue[1]="main"
    if(length(contenu_vue)>1) {
      description_vue[2:length(contenu_vue)]=paste("sup",1:(length(contenu_vue)-1),sep="")
    }
  }
  
  if(!exists("description_vue_analyse_generiq")) {
    description_vue_analyse_generiq=description_vue
  }
  
  n_vue = length(description_vue)
  
  if(!exists("groupe_maximum")) {
    groupe_maximum=list()
    for(i in 1:length(contenu_vue)) {
      groupe_maximum[[i]] = rep(1,length(contenu_vue[[i]]))
    }
  }
  
  if(!exists("repeated")) repeated=FALSE
}

##creation des différents dossiers
cat("Create folders\n")

if(!exists(repertoire)){
  dir.create(repertoire)
}

if(!dir.exists(paste0(repertoire, "/images"))){
  dir.create(paste0(repertoire, "/images"))
}

if(!dir.exists(paste0(repertoire, "/images/boutons"))){
  dir.create(paste0(repertoire, "/images/boutons"))
}

if(!dir.exists(paste0(repertoire, "/analyse"))){
  dir.create(paste0(repertoire, "/analyse"))
}

for(file in c("flechegauche.gif","flechedroite.gif",
              "bouton_simple.gif", "SGD_logo.gif", "traitement.php",
              "design.css")){
  file.copy(from = paste0(common_files_dir, file), 
            to = paste0(repertoire, "/images/boutons/"))
}

file.copy(from = paste0(common_files_dir, c("traitement.php","design.css")), 
          to = repertoire)
file.copy(from = paste0(common_files_dir, c("design.css")), 
          to = paste0(repertoire, "/images"))
file.copy(from = paste0(common_files_dir, c("design.css")), 
          to = paste0(repertoire, "/analyse"))

createIndexHTM(titre_page,repertoire,description_vue,contenu_vue,typeVisu,description_data,decal,texte_description,annot,display_name)

cat("Make count tables\n")

annot_noncount=annot[!(annot$Type %in% noncount),]
nReadsByTranscriptsAnalyse=countReadsBams(data,description_data,annot_noncount,ncore=ncore,normalized=FALSE)
nReadsByTranscriptsNormalAnalyse=countReadsBams(data,description_data,annot_noncount,ncore=ncore,normalized=TRUE)

##creation des fichiers html
S = paste(repertoire,"/analyse/index.htm",sep="")
newHTML(S,paste("Generic analysis for ",titre_page,sep=""))
cat("<p><a href=\"../index.htm\">Return home</a></p>\n",file=S,append=TRUE)
cat("<ul>\n",file=S,append=TRUE)
cat("<li><a href=\"readcount_brut.tab\">transcripts raw read counts</a></li>\n",file=S,append=TRUE)
write.table(cbind(annot_noncount,nReadsByTranscriptsAnalyse),file=paste(repertoire,"/analyse/readcount_brut.tab",sep=""),quote=FALSE,sep="\t",row.names=FALSE)
cat("<li><a href=\"readcount_normal.tab\">transcripts normalized read counts</a></li>\n",file=S,append=TRUE)
write.table(cbind(annot_noncount,nReadsByTranscriptsNormalAnalyse),file=paste(repertoire,"/analyse/readcount_normal.tab",sep=""),quote=FALSE,sep="\t",row.names=FALSE)
nom_fichier_html = "mapping_stats.htm"
S2 = paste(repertoire,"/analyse/",nom_fichier_html,sep="")
newHTML(S2,paste("General mapping stats for ",titre_page,sep=""))
cat("<p><a href=\"index.htm\">Return</a></p>\n",file=S2,append=TRUE)
table_mapping_stats=data.frame(case=names(data[[1]]$flagstat),stringsAsFactors = FALSE)
for(i in 1:length(data)) {
  table_mapping_stats = cbind(table_mapping_stats,data.frame(data[[i]]$flagstat,stringsAsFactors = FALSE))
}
names(table_mapping_stats)[(1:length(data))+1]=description_data
cat("<div class=\"analyse\">\n",file=S2,append=TRUE)
addTable(S2,table_mapping_stats)
cat("</div>\n",file=S2,append=TRUE)
cat("<p><a href=\"mapping_stats.tab\">Download as tab file</a></p>\n",file=S2,append=TRUE)
write.table(table_mapping_stats,file=paste(repertoire,"/analyse/mapping_stats.tab",sep=""),quote=FALSE,sep="\t",row.names=FALSE)
closeHTML(S2)
cat(paste("<li><a href=\"",nom_fichier_html,"\">General mapping stats</a></li>\n",sep=""),file=S,append=TRUE)


##cas des mappings multiples
if(repeated==TRUE) {
  ##creation fichiers html
  nom_fichier_html = "mapping_stats_repeated.htm"
  S2 = paste(repertoire,"/analyse/",nom_fichier_html,sep="")
  newHTML(S2,paste("General mapping stats for ",titre_page," on repeated sequences",sep=""))
  cat("<p><a href=\"index.htm\">Return</a></p>\n",file=S2,append=TRUE)
  cat("<ul>\n",file=S2,append=TRUE)
  for(j in 1:(dim(description_repeated)[1])) {
    nom_fichier_html_ = paste("mapping_stats_repeated_",description_repeated$name[j],".htm",sep="")
    S3 = paste(repertoire,"/analyse/",nom_fichier_html_,sep="")
    newHTML(S3,paste("General mapping stats for ",titre_page," on ",description_repeated$name[j],sep=""))
    cat(paste("<p><a href=\"",nom_fichier_html,"\">Return</a></p>\n",sep=""),file=S3,append=TRUE)
    data_tmp = data_repeated[[j]]
    table_mapping_stats=data.frame(case=names(data[[1]]$flagstat),stringsAsFactors = FALSE)
    for(i in 1:length(data)) {
      if(class(data_tmp[[i]])!="try-error") {
	table_mapping_stats = cbind(table_mapping_stats,data.frame(data_tmp[[i]]$flagstat,stringsAsFactors = FALSE))
      }
      else {
	table_mapping_stats = cbind(table_mapping_stats,rep(NA,length(data[[1]]$flagstat)))
      }
    }
    names(table_mapping_stats)[(1:length(data))+1]=description_data
    cat("<div class=\"analyse\">\n",file=S3,append=TRUE)
    addTable(S3,table_mapping_stats)
    cat("</div>\n",file=S3,append=TRUE)
    nom_fichier_mapping_stats = paste("mapping_stats_",description_repeated$name[j],".tab",sep="")
    cat(paste("<p><a href=\"",nom_fichier_mapping_stats,"\">Download as tab file</a></p>\n"),file=S3,append=TRUE)
    write.table(table_mapping_stats,file=paste(repertoire,"/analyse/",nom_fichier_mapping_stats,sep=""),quote=FALSE,sep="\t",row.names=FALSE)
    closeHTML(S3)
    cat(paste("<li><a href=\"",nom_fichier_html_,"\">General mapping stats for ",description_repeated$name[j],"</a></li>\n",sep=""),file=S2,append=TRUE)
  }
  cat(paste("<li><a href=\"",nom_fichier_html,"\">General mapping stats on repeated sequences</a></li>\n",sep=""),file=S,append=TRUE)
  closeHTML(S2)
}

cat("</ul>\n\n",file=S,append=TRUE)

cat("Make analysis files\n")

nombre_de_transcrits = c()
for(type in unique(annot_noncount$Type)) {
  nombre_de_transcrits=c(nombre_de_transcrits,sum(annot_noncount$Type==type))
}
type_transcrits_tries = unique(annot_noncount$Type)[order(nombre_de_transcrits,decreasing=TRUE)]

taille_figures = c(480,1000)
alias_taille_figures = c("SD","HD")

cat("<div class=\"analyse\">\n",file=S,append=TRUE)
for(vue in 1:length(contenu_vue_analyse_generiq)) {
  contenu_vue_analyse_generiq_ = contenu_vue_analyse_generiq[[vue]]
  n_data_in_vue = length(contenu_vue_analyse_generiq_)
  cat(paste("<h2>",description_vue_analyse_generiq[vue],"</h2>\n",sep=""),file=S,append=TRUE)
  cat("<ul>\n",file=S,append=TRUE)
  nom_fichier_html = paste("mapping_stats_",vue,".htm",sep="")
  S2 = paste(repertoire,"/analyse/",nom_fichier_html,sep="")
  newHTML(S2,paste("Mapping stats ",description_vue_analyse_generiq[vue]," for ",titre_page,sep=""))
  cat("<p><a href=\"index.htm\">Return</a></p>\n",file=S2,append=TRUE)
  table_mapping_stats=data.frame(case=names(data[[1]]$flagstat),stringsAsFactors = FALSE)
  for(i in contenu_vue_analyse_generiq_) {
    table_mapping_stats = cbind(table_mapping_stats,data.frame(data[[i]]$flagstat,stringsAsFactors = FALSE))
  }
  names(table_mapping_stats)[(1:n_data_in_vue)+1]=description_data[contenu_vue_analyse_generiq_]
  cat("<div class=\"analyse\">\n",file=S2,append=TRUE)
  addTable(S2,table_mapping_stats)
  cat("</div>\n",file=S2,append=TRUE)
  nom_fichier_tab = paste("mapping_stats_",vue,".tab",sep="")
  cat(paste("<p><a href=\"",nom_fichier_tab,"\">Download as tab file</a></p>\n",sep=""),file=S2,append=TRUE)
  write.table(table_mapping_stats,file=paste(repertoire,"/analyse/",nom_fichier_tab,sep=""),quote=FALSE,sep="\t",row.names=FALSE)
  closeHTML(S2)
  cat(paste("<li><a href=\"",nom_fichier_html,"\">Mapping stats</a></li>\n",sep=""),file=S,append=TRUE)
  if(repeated==TRUE) {
    nom_fichier_html = paste("mapping_stats_",vue,"_repeated.htm",sep="")
    S2 = paste(repertoire,"/analyse/",nom_fichier_html,sep="")
    newHTML(S2,paste("Mapping stats ",description_vue_analyse_generiq[vue]," for ",titre_page," on repeated sequences",sep=""))
    cat("<p><a href=\"index.htm\">Return</a></p>\n",file=S2,append=TRUE)
    cat("<ul>\n",file=S2,append=TRUE)
    for(j in 1:(dim(description_repeated)[1])) {
      nom_fichier_html_ = paste("mapping_stats_",vue,"_repeated_",description_repeated$name[j],".htm",sep="")
      S3 = paste(repertoire,"/analyse/",nom_fichier_html_,sep="")
      newHTML(S3,paste("Mapping stats ",description_vue_analyse_generiq[vue]," for ",titre_page," on ",description_repeated$name[j],sep=""))
      cat(paste("<p><a href=\"",nom_fichier_html,"\">Return</a></p>\n",sep=""),file=S3,append=TRUE)
      data_tmp = data_repeated[[j]]
      table_mapping_stats=data.frame(case=names(data[[1]]$flagstat),stringsAsFactors = FALSE)
      for(i in contenu_vue_analyse_generiq_) {
	if(class(data_tmp[[i]])!="try-error") {
	  table_mapping_stats = cbind(table_mapping_stats,data.frame(data_tmp[[i]]$flagstat,stringsAsFactors = FALSE))
	}
	else {
	  table_mapping_stats = cbind(table_mapping_stats,rep(NA,length(data[[1]]$flagstat)))
	}
      }
      names(table_mapping_stats)[(1:n_data_in_vue)+1]=description_data[contenu_vue_analyse_generiq_]
      cat("<div class=\"analyse\">\n",file=S3,append=TRUE)
      addTable(S3,table_mapping_stats)
      cat("</div>\n",file=S3,append=TRUE)
      nom_fichier_tab = paste("mapping_stats_",vue,"_",description_repeated$name[j],".tab",sep="")
      cat(paste("<p><a href=\"",nom_fichier_tab,"\">Download as tab file</a></p>\n",sep=""),file=S3,append=TRUE)
      write.table(table_mapping_stats,file=paste(repertoire,"/analyse/",nom_fichier_tab,sep=""),quote=FALSE,sep="\t",row.names=FALSE)
      closeHTML(S3)
      cat(paste("<li><a href=\"",nom_fichier_html_,"\">Mapping stats ",description_vue_analyse_generiq[vue]," for ",description_repeated$name[j],"</a></li>\n",sep=""),file=S2,append=TRUE)
    }
    cat("</ul>\n",file=S2,append=TRUE)
    closeHTML(S2)
    cat(paste("<li><a href=\"",nom_fichier_html,"\">Mapping stats on repeated</a></li>\n",sep=""),file=S,append=TRUE)    
  }
  if(n_data_in_vue*n_data_in_vue<500) {
    nom_fichier_html = paste("scatter_plots_",vue,".htm",sep="")
    S2 = paste(repertoire,"/analyse/",nom_fichier_html,sep="")
    newHTML(S2,paste("Scatter plots ",description_vue_analyse_generiq[vue]," for ",titre_page,sep=""))
    cat("<p><a href=\"index.htm\">Return</a></p>\n",file=S2,append=TRUE)
    cat("<ul>\n",file=S2,append=TRUE)
    top = NA
    bottom = NA
    for(i in 1:n_data_in_vue) {
      top = max(top,nReadsByTranscriptsNormalAnalyse[,2*contenu_vue_analyse_generiq_[i]-1],na.rm=TRUE)
      tmp = nReadsByTranscriptsNormalAnalyse[,2*contenu_vue_analyse_generiq_[i]-1]
      bottom = min(bottom,tmp[tmp>0],na.rm=TRUE)
    }
    cat("<li>All RNAs",file=S2,append=TRUE)
    for(k in 1:length(taille_figures)) {
      taille_figures_ = taille_figures[k]
      magnifier=taille_figures_/480
      chemin_figure = paste("scatter_plot_all_RNAs_",vue,"_",alias_taille_figures[k],".png",sep="")
      cat(paste(" - <a href=\"",chemin_figure,"\">",alias_taille_figures[k],"</a>",sep=""),file=S2,append=TRUE)
      png(paste(repertoire,"/analyse/",chemin_figure,sep=""),taille_figures_,taille_figures_)
      layout(matrix(1:(n_data_in_vue^2),n_data_in_vue,n_data_in_vue))
      for(i in 1:n_data_in_vue) {
	par(lwd=magnifier,mar=c(0.25,0.25,0.25,0.25),xaxt="n",yaxt="n")
	for(j in 1:n_data_in_vue) {
	  if(i!=j) {
	    ind1=contenu_vue_analyse_generiq_[i]
	    ind2=contenu_vue_analyse_generiq_[j]
	    density_cond1 = nReadsByTranscriptsNormalAnalyse[,2*ind1-1]
	    density_cond2 = nReadsByTranscriptsNormalAnalyse[,2*ind2-1]
	    X=log2(density_cond1)
	    Y=log2(density_cond2)
	    X[is.infinite(X)]=NA
	    Y[is.infinite(Y)]=NA
	    bornes = log2(c(bottom,top))
	    plot(X,Y,xlim=bornes,ylim=bornes,main="",xlab="",ylab="",pch=19,cex=0.5)
	    abline(a=0,b=1,col=2)
	  }
	  else {
	    plot(0:2,0:2,col="white",main="",xlab="",ylab="",bty="n")
	    text(1,1,description_data[contenu_vue_analyse_generiq_[i]])
	  }
	}
      }
      dev.off()
    }
    if(length(type_transcrits_tries)>1) {
      for(i_type in 1:length(type_transcrits_tries)) {
	cat(paste("<li>",type_transcrits_tries[i_type],sep=""),file=S2,append=TRUE)
	for(k in 1:length(taille_figures)) {
	  taille_figures_ = taille_figures[k]
	  magnifier=taille_figures_/480
	  chemin_figure = paste("scatter_plot_all_RNAs_",i_type,"_",vue,"_",alias_taille_figures[k],".png",sep="")
	  cat(paste(" - <a href=\"",chemin_figure,"\">",alias_taille_figures[k],"</a>",sep=""),file=S2,append=TRUE)
	  png(paste(repertoire,"/analyse/",chemin_figure,sep=""),taille_figures_,taille_figures_)
	  layout(matrix(1:(n_data_in_vue^2),n_data_in_vue,n_data_in_vue))
	  for(i in 1:n_data_in_vue) {
	    par(lwd=magnifier,mar=c(0.25,0.25,0.25,0.25),xaxt="n",yaxt="n")
	    for(j in 1:n_data_in_vue) {
	      if(i!=j) {
		ind1=contenu_vue_analyse_generiq_[i]
		ind2=contenu_vue_analyse_generiq_[j]
		density_cond1 = nReadsByTranscriptsNormalAnalyse[,2*ind1-1]
		density_cond2 = nReadsByTranscriptsNormalAnalyse[,2*ind2-1]
		X=log2(density_cond1)
		Y=log2(density_cond2)
		X[is.infinite(X)]=NA
		Y[is.infinite(Y)]=NA
		bornes = log2(c(bottom,top))
		plot(X[annot_noncount$Type==type_transcrits_tries[i_type]],Y[annot_noncount$Type==type_transcrits_tries[i_type]],xlim=bornes,ylim=bornes,main="",xlab="",ylab="",pch=19,cex=0.5)
		abline(a=0,b=1,col=2)
	      }
	      else {
		plot(0:2,0:2,col="white",main="",xlab="",ylab="",bty="n")
		text(1,1,description_data[contenu_vue_analyse_generiq_[i]])
	      }
	    }
	  }
	  dev.off()
	}
      }
    }  
    cat(paste("<li><a href=\"",nom_fichier_html,"\">Scatter plots</a></li>",sep=""),file=S,append=TRUE)
  }
  cat("</ul>\n",file=S,append=TRUE)

  cat("<table>\n",file=S,append=TRUE)
  for(i in 1:n_data_in_vue) {
    cat("<tr>\n",file=S,append=TRUE)    
    for(j in 1:n_data_in_vue) {
      if(i!=j) {
	ind1=contenu_vue_analyse_generiq_[i]
	ind2=contenu_vue_analyse_generiq_[j]
	lab_cond1 =description_data[ind1]
	lab_cond2 =description_data[ind2]
	label_abscisses = paste(lab_cond1," log2 densities",sep="")
	label_ordonnees = paste(lab_cond2," log2 densities",sep="")
	density_cond1 = nReadsByTranscriptsNormalAnalyse[,2*ind1-1]
	density_cond2 = nReadsByTranscriptsNormalAnalyse[,2*ind2-1]
	nReads_cond1 = nReadsByTranscriptsNormalAnalyse[,2*ind1]
	nReads_cond2 = nReadsByTranscriptsNormalAnalyse[,2*ind2]
	cond1_over_cond2 = (nReads_cond1+1)/(nReads_cond2+1)
	cond2_over_cond1 = (nReads_cond2+1)/(nReads_cond1+1)
	log_ratio=log2(cond2_over_cond1)
	X=log2(density_cond1)
	Y=log2(density_cond2)
	X[is.infinite(X)]=NA
	Y[is.infinite(Y)]=NA
	correlation=cor(X,Y,use="complete.obs")
	titre_figure = paste(description_data[ind2]," vs ",description_data[ind1],sep="")
        nom_fichier_html = paste("analyse_",ind1,"_",ind2,".htm",sep="")
	S2 = paste(repertoire,"/analyse/",nom_fichier_html,sep="")
	newHTML(S2,paste("Generic analysis for ",titre_figure,sep=""))
        cat("<p><a href=\"index.htm\">Return</a></p>\n",file=S2,append=TRUE)
	cat("\n<ul>\n",file=S2,append=TRUE)
        cat("<li>2D density",file=S2,append=TRUE)
	for(k in 1:length(taille_figures)) {
	  taille_figures_ = taille_figures[k]
	  chemin_figure = paste("2D_plot_densites_",ind1,"_",ind2,"_",alias_taille_figures[k],".png",sep="")
	  cat(paste(" - <a href=\"",chemin_figure,"\">",alias_taille_figures[k],"</a>",sep=""),file=S2,append=TRUE)
	  png(paste(repertoire,"/analyse/",chemin_figure,sep=""),taille_figures_,taille_figures_)
	  magnifier=taille_figures_/480
	  par(lwd=magnifier,cex=magnifier)
	  plot(X,Y,main=titre_figure,xlab=label_abscisses,ylab=label_ordonnees)
	  abline(a=0,b=1,col=2)
	  axis(1,labels=FALSE,lwd.ticks=magnifier)
	  axis(2,labels=FALSE,lwd.ticks=magnifier)
	  mtext(paste("Pearson correlation = ",round(correlation,digits=2)))
	  dev.off()
	}
        cat("</li>\n\n",file=S2,append=TRUE)
	if(length(type_transcrits_tries)>1) {
	  cat("<li>2D density coloured by type",file=S2,append=TRUE)
	  for(k in 1:length(taille_figures)) {
	    nombre_couleur = length(type_transcrits_tries)-1
	    palette=(2:8)[(((1:nombre_couleur)-1)%%7)+1]
	    typepoints=(((1:nombre_couleur)-1)%/%7)+1
	    taille_figures_ = taille_figures[k]
	    chemin_figure = paste("2D_plot_densites_",ind1,"_",ind2,"_",alias_taille_figures[k],"_by_type.png",sep="")
	    cat(paste(" - <a href=\"",chemin_figure,"\">",alias_taille_figures[k],"</a>",sep=""),file=S2,append=TRUE)
	    png(paste(repertoire,"/analyse/",chemin_figure,sep=""),taille_figures_,taille_figures_)
	    magnifier=taille_figures_/480
	    par(lwd=magnifier,cex=magnifier)
	    plot(X,Y,main=titre_figure,xlab=label_abscisses,ylab=label_ordonnees)
	    for(i_type in 2:length(type_transcrits_tries)) {
	      points(X[annot_noncount$Type==type_transcrits_tries[i_type]],Y[annot_noncount$Type==type_transcrits_tries[i_type]],col=palette[i_type-1],pch=typepoints[i_type-1])
	    }
	    legend(x="topleft",legend=type_transcrits_tries,col=c(1,palette),pch=c(1,typepoints),bty="n")
	    abline(a=0,b=1,col=2)
	    axis(1,labels=FALSE,lwd.ticks=magnifier)
	    axis(2,labels=FALSE,lwd.ticks=magnifier)
	    mtext(paste("Pearson correlation = ",round(correlation,digits=2)))
	    dev.off()
	  }
	  cat("</li>\n\n",file=S2,append=TRUE)
	  borne_bas_x=min(X,na.rm=TRUE)
	  borne_haut_x=max(X,na.rm=TRUE)
	  borne_bas_y=min(Y,na.rm=TRUE)
	  borne_haut_y=max(Y,na.rm=TRUE)
	  for(i_type in 2:length(type_transcrits_tries)) {
	    cat(paste("<li>2D density ",type_transcrits_tries[1]," & ",type_transcrits_tries[i_type],sep=""),file=S2,append=TRUE)
	    for(k in 1:length(taille_figures)) {
	      taille_figures_ = taille_figures[k]
	      chemin_figure = paste("2D_plot_densites_",ind1,"_",ind2,"_type_",i_type,"_",alias_taille_figures[k],"_by_type.png",sep="")
	      cat(paste(" - <a href=\"",chemin_figure,"\">",alias_taille_figures[k],"</a>",sep=""),file=S2,append=TRUE)
	      png(paste(repertoire,"/analyse/",chemin_figure,sep=""),taille_figures_,taille_figures_)
	      magnifier=taille_figures_/480
	      par(lwd=magnifier,cex=magnifier)
	      plot(X[annot_noncount$Type==type_transcrits_tries[1]],Y[annot_noncount$Type==type_transcrits_tries[1]],main=titre_figure,xlab=label_abscisses,ylab=label_ordonnees,col=8,xlim=c(borne_bas_x,borne_haut_x),ylim=c(borne_bas_y,borne_haut_y))
	      points(X[annot_noncount$Type==type_transcrits_tries[i_type]],Y[annot_noncount$Type==type_transcrits_tries[i_type]],col=2)
	      legend(x="topleft",legend=c(type_transcrits_tries[1],type_transcrits_tries[i_type]),col=c(8,2),pch=1,bty="n")
	      abline(a=0,b=1,col=2)
	      axis(1,labels=FALSE,lwd.ticks=magnifier)
	      axis(2,labels=FALSE,lwd.ticks=magnifier)
	      dev.off()
	    }
	    cat("</li>\n\n",file=S2,append=TRUE)
	  }
	  breaks_hist = seq(from=min(log_ratio),to=max(log_ratio),length.out=31)
	  hist_list = list()
	  top_hist = 0
	  density_list = list()
	  top_density = 0
	  sd_list = numeric()
	  for(i_type in 1:length(type_transcrits_tries)) {
	    hist_list[[i_type]]=hist(log_ratio[annot_noncount$Type==type_transcrits_tries[i_type]],plot=FALSE,breaks=breaks_hist)
	    top_hist = max(top_hist,hist_list[[i_type]]$density)
	    density_list[[i_type]]=density(log_ratio[annot_noncount$Type==type_transcrits_tries[i_type]],from=min(log_ratio),to=max(log_ratio))
	    top_density=max(top_density,density_list[[i_type]]$y)
	    sd_list[i_type]=sd(log_ratio[annot_noncount$Type==type_transcrits_tries[i_type]])
	  }
	  cat("<li>Histogram number read ratio by type",file=S2,append=TRUE)
	  for(k in 1:length(taille_figures)) {
	    nombre_couleur = length(type_transcrits_tries)-1
	    palette=c(1,(2:8)[(((1:nombre_couleur)-1)%%7)+1])
	    typeligne=c(1,(((1:nombre_couleur)-1)%/%7)+1)
	    taille_figures_ = taille_figures[k]
	    chemin_figure = paste("hist_ratio_",ind1,"_",ind2,"_",alias_taille_figures[k],".png",sep="")
	    cat(paste(" - <a href=\"",chemin_figure,"\">",alias_taille_figures[k],"</a>",sep=""),file=S2,append=TRUE)
	    png(paste(repertoire,"/analyse/",chemin_figure,sep=""),taille_figures_,taille_figures_)
	    magnifier=taille_figures_/480
	    par(lwd=magnifier,cex=magnifier)
	    plot(c(min(log_ratio),max(log_ratio)),c(0,top_hist),main=paste("number read ratio ",lab_cond2,"/",lab_cond1,sep=""),xlab="log2 ratio",ylab="density",col="white")
	    for(i_type in 1:length(type_transcrits_tries)) {
	      lines(hist_list[[i_type]]$mids,hist_list[[i_type]]$density,col=palette[i_type],lty=typeligne[i_type])
	    }
	    axis(1,labels=FALSE,lwd.ticks=magnifier)
	    axis(2,labels=FALSE,lwd.ticks=magnifier)
	    legend(x="topleft",legend=paste(type_transcrits_tries," (sd=",round(sd_list,digits=2),")",sep=""),col=palette,lty=typeligne)
	    dev.off()
	  }
	  cat("</li>\n\n",file=S2,append=TRUE)
	  cat("<li>Density plot read ratio by type",file=S2,append=TRUE)
	  for(k in 1:length(taille_figures)) {
	    nombre_couleur = length(type_transcrits_tries)-1
	    palette=c(1,(2:8)[(((1:nombre_couleur)-1)%%7)+1])
	    typeligne=c(1,(((1:nombre_couleur)-1)%/%7)+1)
	    taille_figures_ = taille_figures[k]
	    chemin_figure = paste("density_ratio_",ind1,"_",ind2,"_",alias_taille_figures[k],".png",sep="")
	    cat(paste(" - <a href=\"",chemin_figure,"\">",alias_taille_figures[k],"</a>",sep=""),file=S2,append=TRUE)
	    png(paste(repertoire,"/analyse/",chemin_figure,sep=""),taille_figures_,taille_figures_)
	    magnifier=taille_figures_/480
	    par(lwd=magnifier,cex=magnifier)
	    plot(c(min(log_ratio),max(log_ratio)),c(0,top_density),main=paste("number read ratio ",lab_cond2,"/",lab_cond1,sep=""),xlab="log2 ratio",ylab="density",col="white")
	    for(i_type in 1:length(type_transcrits_tries)) {
	      lines(density_list[[i_type]]$x,density_list[[i_type]]$y,col=palette[i_type],lty=typeligne[i_type])
	    }
	    axis(1,labels=FALSE,lwd.ticks=magnifier)
	    axis(2,labels=FALSE,lwd.ticks=magnifier)
	    legend(x="topleft",legend=paste(type_transcrits_tries," (sd=",round(sd_list,digits=2),")",sep=""),col=palette,lty=typeligne,bty="n")
	    dev.off()
	  }
	  cat("</li>\n\n",file=S2,append=TRUE)
	}
	tab_with_FC = data.frame(annot_noncount[,1:7],round(nReads_cond1,digits=1),round(nReads_cond2,digits=1),round(cond1_over_cond2,digits=1),round(cond2_over_cond1,digits=1))
	names(tab_with_FC)[8:11]=c(lab_cond1,lab_cond2,paste(lab_cond1,"over",lab_cond2),paste(lab_cond2,"over",lab_cond1))
	FC_moyen_cat = numeric(length(type_transcrits_tries))
	for(i_type in 1:length(type_transcrits_tries)) {
	  FC_moyen_cat[i_type]=mean(cond2_over_cond1[annot_noncount$Type==type_transcrits_tries[i_type]])
	}
	tab_FC_moyen_cat=data.frame(type=type_transcrits_tries,mean_FC=FC_moyen_cat,stringsAsFactors=FALSE)
	names(tab_FC_moyen_cat)[2]=paste("mean",lab_cond2,"over",lab_cond1)
	nom_FC_moyen_file = paste("FC_moyen_",gsub(" ","_",lab_cond1),"_",gsub(" ","_",lab_cond2),".tab",sep="")
	nom_FC_html_file = paste("FC_moyen_",ind1,"_",ind2,".htm",sep="")
	S3 = paste(repertoire,"/analyse/",nom_FC_html_file,sep="")
	cat("<li><a href=\"",nom_FC_html_file,"\">Mean fold-change</a></li>",file=S2,append=TRUE)
	newHTML(S3,paste("Mean fold-change for ",lab_cond2," vs ",lab_cond1,sep=""))
	cat(paste("<p><a href=\"",nom_fichier_html,"\">Return</a></p>\n",sep=""),file=S3,append=TRUE)
	cat("<div class=\"analyse\">\n",file=S3,append=TRUE)
	addTable(S3,tab_FC_moyen_cat)
	cat("</div>\n",file=S3,append=TRUE)
	cat(paste("<p><a href=\"",nom_FC_moyen_file,"\">Download as tab file</a></p>\n",sep=""),file=S3,append=TRUE)
	write.table(tab_FC_moyen_cat,file=paste(repertoire,"/analyse/",nom_FC_moyen_file,sep=""),quote=FALSE,sep="\t",row.names=FALSE)
	closeHTML(S3)
	max_log_ratio = as.integer(max(log_ratio))
	min_log_ratio = as.integer(min(log_ratio))
	list_FC = min_log_ratio:max_log_ratio
	nombre_transcrit_par_type_FC=as.data.frame(matrix(0,length(list_FC),length(type_transcrits_tries)))
	display_FC = numeric(length(list_FC))
	for(i_FC in 1:length(list_FC)) {
	  if(list_FC[i_FC]>0) {
	    display_FC[i_FC]=paste("≥ ",2^list_FC[i_FC],sep="")
	    for(i_type in 1:length(type_transcrits_tries)) {
	      nombre_transcrit_par_type_FC[i_FC,i_type]=sum(log_ratio[annot_noncount$Type==type_transcrits_tries[i_type]]>=list_FC[i_FC])
	    }
	  }
	  if(list_FC[i_FC]<0) {
	    display_FC[i_FC]=paste("≤ 1/",2^(-list_FC[i_FC]),sep="")
	    for(i_type in 1:length(type_transcrits_tries)) {
	      nombre_transcrit_par_type_FC[i_FC,i_type]=sum(log_ratio[annot_noncount$Type==type_transcrits_tries[i_type]]<=list_FC[i_FC])
	    }
	  }
	  if(list_FC[i_FC]==0) {
	    display_FC[i_FC]=">1/2 & <2"
	    for(i_type in 1:length(type_transcrits_tries)) {
	      log_ratio_type=log_ratio[annot_noncount$Type==type_transcrits_tries[i_type]]
	      nombre_transcrit_par_type_FC[i_FC,i_type]=sum(log_ratio_type>-1&log_ratio_type<1)
	    }
	  }
	}
	n_total_transcript_by_type = numeric(length(type_transcrits_tries))
	for(i_type in 1:length(type_transcrits_tries)) {
	  n_total_transcript_by_type[i_type]=sum(annot_noncount$Type==type_transcrits_tries[i_type])
	}
	tab_FC_type=data.frame(c(display_FC,"total"),rbind(nombre_transcrit_par_type_FC,n_total_transcript_by_type),stringsAsFactors=FALSE)
	names(tab_FC_type)=c(paste(lab_cond2,"over",lab_cond1),type_transcrits_tries)
	nom_FC_type_file = paste("FC_type_",gsub(" ","_",lab_cond1),"_",gsub(" ","_",lab_cond2),".tab",sep="")
	nom_FC_type_html_file = paste("FC_type_",ind1,"_",ind2,".htm",sep="")
	S3 = paste(repertoire,"/analyse/",nom_FC_type_html_file,sep="")
	cat("<li><a href=\"",nom_FC_type_html_file,"\">Fold-change thresholds</a></li>",file=S2,append=TRUE)
	newHTML(S3,paste("Fold-change thresholds for ",lab_cond2," vs ",lab_cond1,sep=""))
	cat(paste("<p><a href=\"",nom_fichier_html,"\">Return</a></p>\n",sep=""),file=S3,append=TRUE)
	cat("<div class=\"analyse\">\n",file=S3,append=TRUE)
	addTable(S3,tab_FC_type)
	cat("</div>\n",file=S3,append=TRUE)
	cat(paste("<p><a href=\"",nom_FC_type_file,"\">Download as tab file</a></p>\n",sep=""),file=S3,append=TRUE)
	write.table(tab_FC_type,file=paste(repertoire,"/analyse/",nom_FC_type_file,sep=""),quote=FALSE,sep="\t",row.names=FALSE)
	closeHTML(S3)
	up_tab = tab_with_FC[order(cond2_over_cond1,decreasing=TRUE)[1:300],]
	down_tab = tab_with_FC[order(cond1_over_cond2,decreasing=TRUE)[1:300],]
	nom_up_tab_file = paste("300_up_",gsub(" ","_",lab_cond1),"_",gsub(" ","_",lab_cond2),".tab",sep="")
	nom_down_tab_file = paste("300_down_",gsub(" ","_",lab_cond1),"_",gsub(" ","_",lab_cond2),".tab",sep="")
	cat("<li><a href=\"",nom_up_tab_file,"\">300 UP</a></li>",file=S2,append=TRUE)
	write.table(up_tab,paste(repertoire,"/analyse/",nom_up_tab_file,sep=""),quote=FALSE,sep="\t",row.names=FALSE)
	cat("<li><a href=\"",nom_down_tab_file,"\">300 DOWN</a></li>",file=S2,append=TRUE)
	write.table(down_tab,paste(repertoire,"/analyse/",nom_down_tab_file,sep=""),quote=FALSE,sep="\t",row.names=FALSE)
	nom_up_down_bytype_html_file = paste("up_down_bytype_",ind1,"_",ind2,".htm",sep="")
	S3 = paste(repertoire,"/analyse/",nom_up_down_bytype_html_file,sep="")
	cat("<li><a href=\"",nom_up_down_bytype_html_file,"\">up & down by type</a></li>",file=S2,append=TRUE)
	newHTML(S3,paste("up & down by transcript type ",lab_cond2," vs ",lab_cond1,sep=""))
	cat(paste("<p><a href=\"",nom_fichier_html,"\">Return</a></p>\n",sep=""),file=S3,append=TRUE)
	cat("<ul>",file=S3,append=TRUE)
	for(i_type in 1:length(type_transcrits_tries)) {
	  tab_with_FC_bytype = tab_with_FC[tab_with_FC$Type==type_transcrits_tries[i_type],]
	  cond2_over_cond1_bytype = cond2_over_cond1[tab_with_FC$Type==type_transcrits_tries[i_type]]
	  cond1_over_cond2_bytype = cond1_over_cond2[tab_with_FC$Type==type_transcrits_tries[i_type]]
	  n_transcrit_bytype = sum(tab_with_FC$Type==type_transcrits_tries[i_type])
	  up_tab_bytype = tab_with_FC_bytype[order(cond2_over_cond1_bytype,decreasing=TRUE)[1:min(300,n_transcrit_bytype)],]
	  down_tab_bytype = tab_with_FC_bytype[order(cond1_over_cond2_bytype,decreasing=TRUE)[1:min(300,n_transcrit_bytype)],]
	  nom_up_tab_file_bytype = paste("300_up_bytype_",type_transcrits_tries[i_type],"_",gsub(" ","_",lab_cond1),"_",gsub(" ","_",lab_cond2),".tab",sep="")
	  nom_down_tab_file_bytype = paste("300_down_bytype_",type_transcrits_tries[i_type],"_",gsub(" ","_",lab_cond1),"_",gsub(" ","_",lab_cond2),".tab",sep="")
	  cat("<li><a href=\"",nom_up_tab_file_bytype,"\">300 UP ",type_transcrits_tries[i_type],"</a></li>",file=S3,append=TRUE)
	  write.table(up_tab_bytype,paste(repertoire,"/analyse/",nom_up_tab_file_bytype,sep=""),quote=FALSE,sep="\t",row.names=FALSE)
	  cat("<li><a href=\"",nom_down_tab_file_bytype,"\">300 DOWN ",type_transcrits_tries[i_type],"</a></li>",file=S3,append=TRUE)
	  write.table(down_tab_bytype,paste(repertoire,"/analyse/",nom_down_tab_file_bytype,sep=""),quote=FALSE,sep="\t",row.names=FALSE)
	}
	cat("</ul>",file=S3,append=TRUE)
        closeHTML(S3)
	cat("</ul>",file=S2,append=TRUE)
        closeHTML(S2)
	cat(paste("<td><a href=\"",nom_fichier_html,"\">",titre_figure,"</a></td>",sep=""),file=S,append=TRUE)
      }
      else {
	cat(paste("<td><h4>",description_data[contenu_vue_analyse_generiq_[[i]]],"</h4></td>\n",sep=""),file=S,append=TRUE)
      }
    }
    cat("</tr>\n",file=S,append=TRUE)    
  }
  cat("</table>\n",file=S,append=TRUE)
}
cat("</div>\n",file=S,append=TRUE)
closeHTML(S)

# stop("NI !")

png(paste(repertoire,"/images/legend.png",sep=""),240,200)
prev=par(no.readonly = TRUE)
par(mar=c(0, 0, 2, 0)+0.1)
plotLegend(style,100,100)
par(prev)
dev.off()

calc_start_window = function(start,stop,taille_fenetre,decal) {
  start_num = (ceiling(((start+stop-taille_fenetre+decal/2)/2)/decal)-1)*decal
  start_num[start_num<=0]=1
  return(format(start_num,scientific=FALSE,trim=TRUE))
}

{
if(exists("flags")) {
  annot_flags=rbind(annot,flags)
}
else {
  annot_flags=annot
}
}

annot_flags=annot_flags[sort(annot_flags$Start,index.return=TRUE)$ix,]
Start = calc_start_window(annot_flags$Start,annot_flags$Stop,taille_fenetre,decal)
Chr = annot_flags$Chr
debut_transcrit=annot_flags$Start
fin_transcrit=annot_flags$Stop
brin_transcrit=rep(2,length(annot_flags$Strand))
brin_transcrit[annot_flags$Strand=="-"]=1
brin_transcrit[annot_flags$Strand=="+"]=0
for(i in 1:n_vue) {
  liens_ = paste("<a href=\"",Chr,"_",Start,"_",i,"_note_",debut_transcrit,"_",fin_transcrit,"_",brin_transcrit,".htm\"><img src=\"boutons/bouton_simple.gif\"/></a>",sep="")
  if(i==1) {
    liens = data.frame(liens_,stringsAsFactors = FALSE)
  }
  else {
    liens = data.frame(liens,liens_,stringsAsFactors = FALSE)
  }
}
infoName=unique(c("ID","gene","Name","Alias",display_name))
N = dim(annot_flags)[1]
gene_=numeric(N)
for(i in 1:N) {
  gene__=unique(unlist(strsplit(as.character(annot_flags[i,names(annot_flags) %in% infoName]),split=",")))
  gene__=gene__[gene__!="."]
  if(length(gene__)==0) gene__="."
  gene_[i]=paste(gene__,collapse="<br/>")
}

SGD_=paste("<a href=\"https://www.yeastgenome.org/locus/",annot_flags$ID,"\"><img src=\"boutons/SGD_logo.gif\"/></a>",sep="")
for(type in nonSGD) {
  SGD_[annot_flags$Type==type]=""
}
names(liens)=description_vue
# disp_annot=data.frame(annot_flags[,c(6,2)],gene_,annot_flags[,3:5],SGD_,liens,stringsAsFactors = FALSE)[!(annot_flags$Type %in% noncount),]
disp_annot=data.frame(annot_flags[,c(6,2)],gene_,annot_flags[,3:5],SGD_,liens,stringsAsFactors = FALSE)

names(disp_annot)[3]="gene"
names(disp_annot)[7]="SGD"

annot_flags_for_count = annot_flags[!(annot_flags$Type %in% noncount),]
disp_annot_for_count = disp_annot[!(annot_flags$Type %in% noncount),]

nReadsByTranscriptsNormal=countReadsBams(data,description_data,annot_flags_for_count,ncore=ncore,normalized=TRUE)
densitiesGenesNormal=format(nReadsByTranscriptsNormal[,(2*1:length(description_data))-1],scientific=TRUE,digits=3)
N=dim(densitiesGenesNormal)[2]
if(is.null(N)) {
  densitiesGenesNormal = array(densitiesGenesNormal,dim=c(1,length(densitiesGenesNormal)))
  N=dim(densitiesGenesNormal)[2]
}
for(i in 1:N) {
  densitiesGenesNormal[,i]=gsub("+","",densitiesGenesNormal[,i],fixed=TRUE)
  densitiesGenesNormal[,i]=gsub("e0","e",densitiesGenesNormal[,i],fixed=TRUE)
  densitiesGenesNormal[,i]=gsub("e-0","e-",densitiesGenesNormal[,i],fixed=TRUE)
  densitiesGenesNormal[,i]=paste(gsub("e"," × 10<sup>",densitiesGenesNormal[,i]),"</sup>",sep="")
  densitiesGenesNormal[,i]=gsub(" × 10<sup>0</sup>","",densitiesGenesNormal[,i],fixed=TRUE)
}

annot_flags_for_form = annot_flags[!(annot_flags$Type %in% nonform),]
disp_annot_for_form = disp_annot[!(annot_flags$Type %in% nonform),]

infoName=unique(c("ID","gene","Name","Alias",display_name))
list_names=list()
for(infoName_ in infoName) {
  if(infoName_ %in% names(annot_flags_for_form)) {
    list_names[[length(list_names)+1]]=strsplit(annot_flags_for_form[,infoName_],split=",",fixed=TRUE)
  }
}

all_pos=numeric()
all_names=numeric()

nombre_gene=dim(annot_flags_for_form)[1]
for(i in 1:nombre_gene) {
  names_=c()
  for(j in 1:length(list_names)) {
    names_ = c(names_,unlist(list_names[[j]][i]))
  }
  names_=unique(names_)
  all_pos = c(all_pos,rep(i,length(names_)))
  all_names = c(all_names,unlist(names_))
}

all_names=toupper(all_names)

makePHPredirect = function(name) {
  pos_ = unique(all_pos[all_names==name])
  if(length(pos_)==1) {
    for(j in 1:n_vue) {
      chr_ = annot_flags_for_form$Chr[pos_]
      start = calc_start_window(annot_flags_for_form$Start[pos_],annot_flags_for_form$Stop[pos_],taille_fenetre,decal)
      brin_=2
      if(annot_flags_for_form$Strand[pos_]=="+") brin_=0
      if(annot_flags_for_form$Strand[pos_]=="-") brin_=1      
      cat(paste("<?php \nheader('Location: ",chr_,"_",start,"_",j,"_note_",annot_flags_for_form$Start[pos_],"_",annot_flags_for_form$Stop[pos_],"_",brin_,".htm');\n ?>",sep=""),file=paste(repertoire,"/images/",name,"_",j,".php",sep=""))
    }
  }
  else {
    nom_fichier=paste(name,"_redirection.htm",sep="")
    for(j in 1:3) {
      cat(paste("<?php \nheader('Location: ",nom_fichier,"');\n ?>",sep=""),file=paste(repertoire,"/images/",name,"_",j,".php",sep=""))
    }
    S=paste(repertoire,"/images/",nom_fichier,sep="")
    newHTML(S,paste(name," redirection",sep=""))
    cat(paste("<p>",name," corresponds to ",length(pos_)," known transcripts.</p>",sep=""),file=S,append=TRUE)
    addTable(S,disp_annot_for_form[pos_,])
    closeHTML(S)
  }
}

cat("Make php redirection files\n")

lapply(unique(all_names[all_names!="."]),makePHPredirect)

convertData = function(chr) {
  forward=list()
  reverse=list()
  chrName=noms_chromosome[chr]
  for(i in 1:length(data)){
    forward_=numeric(tailles_chromosome[chr])
    n_reads = length(data[[i]]$F[[chrName]][[1]])
    if(n_reads>0) {
      for(k in 1:n_reads) {
	debut_read = max(1,data[[i]]$F[[chrName]][[1]][k])
	fin_read = min(data[[i]]$F[[chrName]][[2]][k],tailles_chromosome[chr])
	forward_[debut_read:fin_read]=forward_[debut_read:fin_read]+data[[i]]$F[[chrName]][[4]][k]
      }
    }
    n_junctions = length(data[[i]]$junctions_F[[chrName]][[1]])
    if(n_junctions>0) {
      for(k in 1:n_junctions) {
	debut_junction = max(1,data[[i]]$junctions_F[[chrName]][[1]][k])
	fin_junction = min(data[[i]]$junctions_F[[chrName]][[2]][k],tailles_chromosome[chr])
	forward_[debut_junction:fin_junction]=forward_[debut_junction:fin_junction]-data[[i]]$junctions_F[[chrName]][[4]][k]
      }
    }

    reverse_=numeric(tailles_chromosome[chr])
    n_reads = length(data[[i]]$R[[chrName]][[1]])
    if(n_reads>0) {
      for(k in 1:n_reads) {
	debut_read = max(1,data[[i]]$R[[chrName]][[1]][k])
	fin_read = min(data[[i]]$R[[chrName]][[2]][k],tailles_chromosome[chr])
	reverse_[debut_read:fin_read]=reverse_[debut_read:fin_read]+data[[i]]$R[[chrName]][[4]][k]
      }
    }
    n_junctions = length(data[[i]]$junctions_R[[chrName]][[1]])
    if(n_junctions>0) {
      for(k in 1:n_junctions) {
	debut_junction = max(1,data[[i]]$junctions_R[[chrName]][[1]][k])
	fin_junction = min(data[[i]]$junctions_R[[chrName]][[2]][k],tailles_chromosome[chr])
	reverse_[debut_junction:fin_junction]=reverse_[debut_junction:fin_junction]-data[[i]]$junctions_R[[chrName]][[4]][k]
      }
    }

    forward_[forward_<0]=0
    reverse_[reverse_<0]=0

    if(data[[i]]$stranded) {
      forward[[i]]=forward_
      reverse[[i]]=reverse_
    }
    else {
      forward[[i]]=forward_+reverse_
      reverse[[i]]=rep(0,length(forward[[i]]))
    }
  }
  chr_ = list()
  chr_$F = forward
  chr_$R = reverse
  chr_$name = noms_chromosome[chr]
  return(chr_)
}

chrProcess = function(i_chr,ncore=1) {
  chr=extractSignal(data,noms_chromosome[i_chr],normalized_=TRUE)
  taille_chromosome=length(chr$F[[1]])
  n_windows=max(ceiling(taille_chromosome/decal)-ceiling(taille_fenetre/decal)+1,1)
  start_windows=(0:n_windows)*decal
  stop_windows=start_windows+taille_fenetre
  start_windows[1]=1
  disp_start_windows=format(start_windows,scientific=FALSE,trim=TRUE)
  disp_stop_windows=format(stop_windows,scientific=FALSE,trim=TRUE)
  
  if(!(isVisuStranded)) {
    for(i in 1:length(chr$F)){
      chr$F[[i]]=chr$F[[i]]+chr$R[[i]]
      chr$R[[i]]=double(length(chr$R[[i]]))
    }    
  }
  
  ##permet de gérer taille fenêtre (on affiche 15000nt mais chromo pas forcément multiple de 15000)
  if(isVisuLog) {
    forward=list()
    reverse=list()
    for(i in 1:length(chr$F)){
      forward_ = double(stop_windows[n_windows])
      tmp=log2(1+chr$F[[i]])
      forward_[1:length(tmp)]=tmp
      forward[[i]]=forward_
      reverse_ = double(stop_windows[n_windows])
      tmp=log2(1+chr$R[[i]])
      reverse_[1:length(tmp)]=tmp
      reverse[[i]]=reverse_    
    }
  }else {
    forward=list()
    reverse=list()
    for(i in 1:length(chr$F)){
      forward_ = double(stop_windows[n_windows])
      tmp=chr$F[[i]]
      forward_[1:length(tmp)]=tmp
      forward[[i]]=forward_
      reverse_ = double(stop_windows[n_windows])
      tmp=chr$R[[i]]
      reverse_[1:length(tmp)]=tmp
      reverse[[i]]=reverse_    
    }  
  }
  
  
  if(exists("is_there_any_chipseq")) {
    if(is_there_any_chipseq) {
      for(i in 1:length(index_chipseq)) {
	i_chip = index_chipseq[i]
	i_cont = index_control[i]
	tmp = forward[[i_chip]] - forward[[i_cont]]
	tmp[tmp<=0] = 0
	forward[[i_chip]] = tmp
	tmp = reverse[[i_chip]] - reverse[[i_cont]]
	tmp[tmp<=0] = 0
	reverse[[i_chip]] = tmp
      }
    }
  }
  
  ## Vue chromosome
  debut_vue = 1
  fin_vue = taille_chromosome
  disp_debut_vue = format(debut_vue,scientific=FALSE)
  disp_fin_vue = format(fin_vue,scientific=FALSE)
  doVueWholeChr = function(vue) {
    cat(paste("<?php \nheader('Location: ",noms_chromosome[i_chr],"_",vue,".htm');\n ?>",sep=""),file=paste(repertoire,"/images/",toupper(noms_chromosome[i_chr]),"_",vue,".php",sep=""))
    S=paste(repertoire,"/images/",noms_chromosome[i_chr],"_",vue,".htm",sep="")
    newHTML(S,paste(noms_chromosome[i_chr],":",disp_debut_vue,"..",disp_fin_vue," - ",description_vue[vue],sep=""))
    addNavigatorHautRepeated(S,noms_chromosome[i_chr],vue,description_vue,decal)
    contenu_vue_=contenu_vue[[vue]]
    n_element_vue=length(contenu_vue_)
    max_forward = numeric(n_element_vue)
    max_reverse = numeric(n_element_vue)
    for(element in 1:n_element_vue) {
      i_data = contenu_vue_[element]
      max_forward[element]=max(forward[[i_data]][debut_vue:fin_vue])
      max_reverse[element]=max(reverse[[i_data]][debut_vue:fin_vue])
    }
    groupe_maximum_=groupe_maximum[[vue]]
    for(i_max in unique(groupe_maximum_)) {
      max_forward[groupe_maximum_==i_max]=max(max_forward[groupe_maximum_==i_max])
      max_reverse[groupe_maximum_==i_max]=max(max_reverse[groupe_maximum_==i_max])
    }
    nom_fichier_image=paste(noms_chromosome[i_chr],"_",vue,".png",sep="")
    cat(paste("<p>\n<img src=\"",nom_fichier_image,"\" />\n</p>\n",sep=""),file=S,append=TRUE)
    chemin_fichier_image = paste(repertoire,"/images/",nom_fichier_image,sep="")
    plotVisu(file=chemin_fichier_image,typeVisu=typeVisu[vue],listForward=forward,listReverse=reverse,which=contenu_vue_,
    debut_vue=debut_vue,fin_vue=fin_vue,chr=noms_chromosome[i_chr],annot=annot,style=style,log=isVisuLog,stranded=isVisuStranded,
    top=max_forward,bottom=max_reverse,titres=description_data,name_flags="flags",
    classic_plus_color=classic_plus_color,classic_minus_color=classic_minus_color,smoothLength=smoothLength,
    heatmap_max_color=heatmap_max_color,heatmap_min_color=heatmap_min_color,heatmap_palette_method=heatmap_palette_method,
    lines_samples_colors=lines_samples_colors[contenu_vue_],lines_samples_type_line=lines_samples_type_line[contenu_vue_],display_name=display_name,initialize_label_sizes=initialize_label_sizes)
    cat("</div><hr/>",file=S,append=TRUE)

    annot_ = cbind(disp_annot_for_count[(annot_flags_for_count$Chr==noms_chromosome[i_chr])&(annot_flags_for_count$Start<=fin_vue)&(annot_flags_for_count$Stop>=debut_vue),],
                   densitiesGenesNormal[(annot_flags_for_count$Chr==noms_chromosome[i_chr])&(annot_flags_for_count$Start<=fin_vue)&(annot_flags_for_count$Stop>=debut_vue),contenu_vue_])
    colnames(annot_) <- c(colnames(disp_annot_for_count), colnames(densitiesGenesNormal)[contenu_vue_])
    addTable(S,annot_)
    
    closeHTML(S)
    return(0)
  }
  
  if(ncore==1) {
    ctlVueWholeChr=lapply(1:n_vue,doVueWholeChr)
  }
  else {
    ctlVueWholeChr=mclapply(1:n_vue,doVueWholeChr,mc.cores=ncore,mc.preschedule=FALSE)
  }

  
  ## Par fenêtre
  doVueWindows = function(i) {
    debut_vue=start_windows[i]
    fin_vue=stop_windows[i]
    disp_debut_vue=disp_start_windows[i]
    disp_fin_vue=disp_stop_windows[i]
    for(vue in 1:n_vue) {
      S=paste(repertoire,"/images/",noms_chromosome[i_chr],"_",disp_debut_vue,"_",vue,".htm",sep="")
      newHTML(S,paste(noms_chromosome[i_chr],":",disp_debut_vue,"..",disp_fin_vue," - ",description_vue[vue],sep=""))
      addNavigatorHaut(S,noms_chromosome[i_chr],debut_vue,vue,description_vue,decal)

      if(i==1) {
	addNavigatorBas(S,noms_chromosome[i_chr],debut_vue,0,start_windows[i+1],decal,taille_fenetre,vue,TRUE)
      }
      else {
	if(i==n_windows) {
	  addNavigatorBas(S,noms_chromosome[i_chr],debut_vue,start_windows[i-1],0,decal,taille_fenetre,vue,TRUE)
	}
	else{
	  addNavigatorBas(S,noms_chromosome[i_chr],debut_vue,start_windows[i-1],start_windows[i+1],decal,taille_fenetre,vue,TRUE)
	}
      }

      contenu_vue_=contenu_vue[[vue]]
      n_element_vue=length(contenu_vue_)
      max_forward = numeric(n_element_vue)
      max_reverse = numeric(n_element_vue)
      for(element in 1:n_element_vue) {
	i_data = contenu_vue_[element]
	max_forward[element]=max(forward[[i_data]][debut_vue:fin_vue])
	max_reverse[element]=max(reverse[[i_data]][debut_vue:fin_vue])
      }
      groupe_maximum_=groupe_maximum[[vue]]
      for(i_max in unique(groupe_maximum_)) {
	max_forward[groupe_maximum_==i_max]=max(max_forward[groupe_maximum_==i_max])
	max_reverse[groupe_maximum_==i_max]=max(max_reverse[groupe_maximum_==i_max])
      }
      nom_fichier_image=paste(noms_chromosome[i_chr],"_",disp_debut_vue,"_",vue,".png",sep="")
      cat(paste("<p>\n<img src=\"",nom_fichier_image,"\" />\n</p>\n",sep=""),file=S,append=TRUE)
      chemin_fichier_image=paste(repertoire,"/images/",nom_fichier_image,sep="")
      plotVisu(file=chemin_fichier_image,typeVisu=typeVisu[vue],listForward=forward,listReverse=reverse,
      which=contenu_vue_,debut_vue=debut_vue,fin_vu=fin_vue,chr=noms_chromosome[i_chr],annot=annot,style=style,log=isVisuLog,stranded=isVisuStranded,
      top=max_forward,bottom=max_reverse,titres=description_data,name_flags="flags",
      classic_plus_color=classic_plus_color,classic_minus_color=classic_minus_color,smoothLength=smoothLength,
      heatmap_max_color=heatmap_max_color,heatmap_min_color=heatmap_min_color,heatmap_palette_method=heatmap_palette_method,
      lines_samples_colors=lines_samples_colors[contenu_vue_],lines_samples_type_line=lines_samples_type_line[contenu_vue_],display_name=display_name,initialize_label_sizes=initialize_label_sizes)

      if(i==1) {
	addNavigatorBas(S,noms_chromosome[i_chr],debut_vue,0,start_windows[i+1],decal,taille_fenetre,vue)
      }
      else {
	if(i==n_windows) {
	  addNavigatorBas(S,noms_chromosome[i_chr],debut_vue,start_windows[i-1],0,decal,taille_fenetre,vue)
	}
	else{
	  addNavigatorBas(S,noms_chromosome[i_chr],debut_vue,start_windows[i-1],start_windows[i+1],decal,taille_fenetre,vue)
	}
      }

      cat("</div><hr/>",file=S,append=TRUE)

      annot_ = cbind(disp_annot_for_count[(annot_flags_for_count$Chr==noms_chromosome[i_chr])&(annot_flags_for_count$Start<=fin_vue)&(annot_flags_for_count$Stop>=debut_vue),],
                     densitiesGenesNormal[(annot_flags_for_count$Chr==noms_chromosome[i_chr])&(annot_flags_for_count$Start<=fin_vue)&(annot_flags_for_count$Stop>=debut_vue),contenu_vue_])
      colnames(annot_) <- c(colnames(disp_annot_for_count), colnames(densitiesGenesNormal)[contenu_vue_])
      addTable(S,annot_)
      
      closeHTML(S)
    }
    return(0)
  }

  if(ncore==1) {
    ctlVueWindows=lapply(1:n_windows,doVueWindows)
  }
  else {
    ctlVueWindows=mclapply(1:n_windows,doVueWindows,mc.cores=ncore,mc.preschedule=FALSE)
  }

  ## Par transcrit
  annot_chr = annot_flags[annot_flags$Chr==noms_chromosome[i_chr],]
  n_transcrit=dim(annot_chr)[1]
  doVueTranscrit = function(j) {
    debut_transcrit=annot_chr$Start[j]
    fin_transcrit=annot_chr$Stop[j]
    if(annot_chr$Strand[j]=="+") {
      brin_transcrit=0
    }else {
      if(annot_chr$Strand[j]=="-") {
	brin_transcrit=1
      }else {
	brin_transcrit=2
      }
    }
    debut_vue=max(1,(ceiling(((debut_transcrit+fin_transcrit-taille_fenetre+decal/2)/2)/decal)-1)*decal)
    i=(1:n_windows)[start_windows==debut_vue]
    i_correct=TRUE
    if(length(i)==0) {
      i_correct=FALSE
    }else {
      if(is.na(i)) {
	i_correct=FALSE
      }
    }
    if(i_correct==FALSE) {
      i=n_windows
      debut_vue=start_windows[i]
    }
    fin_vue=stop_windows[i]
    disp_debut_vue=disp_start_windows[i]    
    disp_fin_vue=disp_stop_windows[i]
    
    for(vue in 1:n_vue) {
      S=paste(repertoire,"/images/",noms_chromosome[i_chr],"_",disp_debut_vue,"_",vue,"_note_",debut_transcrit,"_",fin_transcrit,"_",brin_transcrit,".htm",sep="")
      newHTML(S,paste(noms_chromosome[i_chr],":",disp_debut_vue,"..",disp_fin_vue," - ",description_vue[vue],sep=""))
      addNavigatorHaut(S,noms_chromosome[i_chr],debut_vue,vue,description_vue,decal)

      if(i==1) {
	addNavigatorBas(S,noms_chromosome[i_chr],debut_vue,0,start_windows[i+1],decal,taille_fenetre,vue,TRUE)
      }else {
	if(i==n_windows) {
	  addNavigatorBas(S,noms_chromosome[i_chr],debut_vue,start_windows[i-1],0,decal,taille_fenetre,vue,TRUE)
	}else{
	  addNavigatorBas(S,noms_chromosome[i_chr],debut_vue,start_windows[i-1],start_windows[i+1],decal,taille_fenetre,vue,TRUE)
	}
      }

      contenu_vue_=contenu_vue[[vue]]
      n_element_vue=length(contenu_vue_)
      max_forward = numeric(n_element_vue)
      max_reverse = numeric(n_element_vue)
      for(element in 1:n_element_vue) {
	i_data = contenu_vue_[element]
	max_forward[element]=max(forward[[i_data]][debut_vue:fin_vue])
	max_reverse[element]=max(reverse[[i_data]][debut_vue:fin_vue])
      }
      groupe_maximum_=groupe_maximum[[vue]]
      for(i_max in unique(groupe_maximum_)) {
	max_forward[groupe_maximum_==i_max]=max(max_forward[groupe_maximum_==i_max])
	max_reverse[groupe_maximum_==i_max]=max(max_reverse[groupe_maximum_==i_max])
      }
      nom_fichier_image=paste(noms_chromosome[i_chr],"_",disp_debut_vue,"_",vue,"_note_",debut_transcrit,"_",fin_transcrit,"_",brin_transcrit,".png",sep="")
      cat(paste("<p>\n<img src=\"",nom_fichier_image,"\" />\n</p>\n",sep=""),file=S,append=TRUE)
      chemin_fichier_image=paste(repertoire,"/images/",nom_fichier_image,sep="")
      plotVisu(file=chemin_fichier_image,typeVisu=typeVisu[vue],listForward=forward,listReverse=reverse,
      which=contenu_vue_,debut_vue=debut_vue,fin_vue=fin_vue,chr=noms_chromosome[i_chr],annot=annot,style=style,
      top=max_forward,bottom=max_reverse,marks=c(debut_transcrit,fin_transcrit),log=isVisuLog,stranded=isVisuStranded,
      strandMarks=brin_transcrit,titres=description_data,name_flags="flags",
      classic_plus_color=classic_plus_color,classic_minus_color=classic_minus_color,smoothLength=smoothLength,
      heatmap_max_color=heatmap_max_color,heatmap_min_color=heatmap_min_color,heatmap_palette_method=heatmap_palette_method,
      lines_samples_colors=lines_samples_colors[contenu_vue_],lines_samples_type_line=lines_samples_type_line[contenu_vue_],display_name=display_name,initialize_label_sizes=initialize_label_sizes)

      if(i==1) {
	addNavigatorBas(S,noms_chromosome[i_chr],debut_vue,0,start_windows[i+1],decal,taille_fenetre,vue)
      }else {
	if(i==n_windows) {
	  addNavigatorBas(S,noms_chromosome[i_chr],debut_vue,start_windows[i-1],0,decal,taille_fenetre,vue)
	}else{
	  addNavigatorBas(S,noms_chromosome[i_chr],debut_vue,start_windows[i-1],start_windows[i+1],decal,taille_fenetre,vue)
	}
      }

      cat("</div><hr/>",file=S,append=TRUE)

      annot_ = cbind(disp_annot_for_count[(annot_flags_for_count$Chr==noms_chromosome[i_chr])&(annot_flags_for_count$Start<=fin_vue)&(annot_flags_for_count$Stop>=debut_vue),],
                     densitiesGenesNormal[(annot_flags_for_count$Chr==noms_chromosome[i_chr])&(annot_flags_for_count$Start<=fin_vue)&(annot_flags_for_count$Stop>=debut_vue),contenu_vue_])
      colnames(annot_) <- c(colnames(disp_annot_for_count), colnames(densitiesGenesNormal)[contenu_vue_])
      addTable(S,annot_)
      
      closeHTML(S)
    }
    return(0)
  }

  if(ncore==1) {
    ctlVueTranscrit=lapply(1:n_transcrit,doVueTranscrit)
  }
  else {
    ctlVueTranscrit=mclapply(1:n_transcrit,doVueTranscrit,mc.cores=ncore,mc.preschedule=FALSE)
  }

  return(list(ctlVueWholeChr,ctlVueWindows,ctlVueTranscrit))
}

initialize_label_sizes=initialize_parking(annot_flags,typeVisu=unique(typeVisu),taille_fenetre=taille_fenetre,display_name=display_name)

cat("Start producing images\n")

ctl=list()
for(i_chr in 1:n_chromosome) {
  
  cat("    chromosome ",noms_chromosome[i_chr],"\n")
  
  ctl[[i_chr]]=try(chrProcess(i_chr,ncore=ncore),TRUE)
}

##mapping multiple
if(repeated==TRUE) {
  repeatProcess = function(i){
    is_ok = TRUE
    for(j in 1:length(description_data)) {
      is_ok = is_ok & (class(data_repeated[[i]][[j]])!="try-error")
    }
    if(is_ok) {
	  ##rajout normalisation
	  index_norm=4
	  if (!data_repeated[[i]][[j]]$norm){
		index_norm=3
	  }
	  ##
      name=index_sequence[i]
      chr=description_repeated$chr[i]
      debut_vue=description_repeated$start[i]
      fin_vue=description_repeated$stop[i]
      longueur_sequence=fin_vue-debut_vue+1

      forward=list()
      reverse=list()
      for(j in 1:length(description_data)) {
		forward_=numeric(longueur_sequence)
		n_reads = length(data_repeated[[i]][[j]]$F[[1]][[1]])
		if(n_reads>0) {
			for(k in 1:n_reads) {
				debut_read = max(1,data_repeated[[i]][[j]]$F[[1]][[1]][k])
				fin_read = min(data_repeated[[i]][[j]]$F[[1]][[2]][k],longueur_sequence)
				forward_[debut_read:fin_read]=forward_[debut_read:fin_read]+data_repeated[[i]][[j]]$F[[1]][[index_norm]][k]
			}	
		}

		##rajout junctions 
		which_junctions=which((data_repeated[[i]][[j]]$junctions_F[[1]][[2]]>=debut_vue) & (data_repeated[[i]][[j]]$junctions_F[[1]][[1]]<=fin_vue))
		n_junctions=length(which_junctions)
		if(n_junctions>0) {
		  for(k in which_junctions) {
			debut_junction=max(1,data_repeated[[i]][[j]]$junctions_F[[1]][[1]][k])
			fin_junction=min(data_repeated[[i]][[j]]$junctions_F[[1]][[2]][k],longueur_sequence)
			forward_[debut_junction:fin_junction]=forward_[debut_junction:fin_junction]-data_repeated[[i]][[j]]$junctions_F[[1]][[index_norm]][k]
		  }
		}
		##
		reverse_=numeric(longueur_sequence)
		n_reads = length(data_repeated[[i]][[j]]$R[[1]][[1]])
		if(n_reads>0) {
		  for(k in 1:n_reads) {
			debut_read = max(1,data_repeated[[i]][[j]]$R[[1]][[1]][k])
			fin_read = min(data_repeated[[i]][[j]]$R[[1]][[2]][k],longueur_sequence)
			reverse_[debut_read:fin_read]=reverse_[debut_read:fin_read]+data_repeated[[i]][[j]]$R[[1]][[index_norm]][k]
		  }
		}
		##rajout junctions 
		which_junctions=which((data_repeated[[i]][[j]]$junctions_R[[1]][[2]]>=debut_vue) & (data_repeated[[i]][[j]]$junctions_R[[1]][[1]]<=fin_vue))
		n_junctions=length(which_junctions)
		if(n_junctions>0) {
		  for(k in which_junctions) {
			debut_junction=max(1,data_repeated[[i]][[j]]$junctions_R[[1]][[1]][k])
			fin_junction=min(data_repeated[[i]][[j]]$junctions_R[[1]][[2]][k],longueur_sequence)
			reverse_[debut_junction:fin_junction]=reverse_[debut_junction:fin_junction]-data_repeated[[i]][[j]]$junctions_R[[1]][[index_norm]][k]
		  }
		}
		##
		if(data_repeated[[i]][[j]]$stranded) {
		  if(isVisuLog) { ## rajout option log
			forward[[j]]=log2(1+forward_)
			reverse[[j]]=log2(1+reverse_)
		  }else{
			forward[[j]]=forward_
			reverse[[j]]=reverse_
		  }
		}
		else {
		  if (isVisuLog){
			forward[[j]]=log2(1+forward_+reverse_)
		   }else{
			forward[[j]]=forward_+reverse_

		   }
		  reverse[[j]]=numeric(longueur_sequence)
		}
      }

      if(exists("is_there_any_chipseq")) {
		if(is_there_any_chipseq) {
		  for(i in 1:length(index_chipseq)) {
			i_chip = index_chipseq[i]
			i_cont = index_control[i]
			tmp = forward[[i_chip]] - forward[[i_cont]]
			tmp[tmp<=0] = 0
			forward[[i_chip]] = tmp
		  }  
		}
      }

      density = numeric()
      
      if(description_repeated$legend[i]=="same") {
	{
	if(data_repeated[[i]][[j]]$stranded) {
	  if(description_repeated$strand[i]=="+") {
	    brin_transcrit=0
	    for(j in 1:length(description_data)) {
	      density[j]=sum(data_repeated[[i]][[j]]$F[[1]][[index_norm]])/longueur_sequence
	    }
	  }
	  else {
	    brin_transcrit=1
	    for(j in 1:length(description_data)) {
	      density[j]=sum(data_repeated[[i]][[j]]$R[[1]][[index_norm]])/longueur_sequence
	    }
	  }
	}
	else {
	  brin_transcrit=0
	  for(j in 1:length(description_data)) {
	    density[j]=(sum(data_repeated[[i]][[j]]$R[[1]][[index_norm]])/longueur_sequence)+(sum(data_repeated[[i]][[j]]$F[[1]][[index_norm]])/longueur_sequence)
	  }
	}
	}
	density_format=format(density,scientific=TRUE,digits=3)
	density_format=gsub("+","",density_format,fixed=TRUE)
	density_format=gsub("e0","e",density_format,fixed=TRUE)
	density_format=gsub("e-0","e-",density_format,fixed=TRUE)
	density_format=paste(gsub("e"," × 10<sup>",density_format),"</sup>",sep="")
	density_format=gsub(" × 10<sup>0</sup>","",density_format,fixed=TRUE)

	quantif = data.frame(description_data,density_format,stringsAsFactors = FALSE)
	names(quantif)=c("Data","Density")
	
	for(vue in 1:n_vue) {
	  cat(paste("<?php \nheader('Location: ",name,"_",vue,".htm');\n ?>",sep=""),file=paste(repertoire,"/images/",toupper(name),"_",vue,".php",sep=""))
	  S=paste(repertoire,"/images/",name,"_",vue,".htm",sep="")
	  newHTML(S,paste(name,"  ",chr,":",debut_vue,"..",fin_vue," - ",description_vue[vue],sep=""))

	  addNavigatorHautRepeated(S,name,vue,description_vue,decal)

	  contenu_vue_=contenu_vue[[vue]]
	  n_element_vue=length(contenu_vue_)
	  max_forward = numeric(n_element_vue)
	  max_reverse = numeric(n_element_vue)
	  for(element in 1:n_element_vue) {
	    i_data = contenu_vue_[element]
	    max_forward[element]=max(forward[[i_data]])
	    max_reverse[element]=max(reverse[[i_data]])
	  }
	  groupe_maximum_=groupe_maximum[[vue]]
	  for(i_max in unique(groupe_maximum_)) {
	    max_forward[groupe_maximum_==i_max]=max(max_forward[groupe_maximum_==i_max])
	    max_reverse[groupe_maximum_==i_max]=max(max_reverse[groupe_maximum_==i_max])
	  }
	  nom_fichier_image=paste(name,"_",vue,".png",sep="")
	  cat(paste("<p>\n<img src=\"",nom_fichier_image,"\" />\n</p>\n",sep=""),file=S,append=TRUE)
	  chemin_fichier_image= paste(repertoire,"/images/",nom_fichier_image,sep="")
	  plotVisu(chemin_fichier_image,typeVisu=typeVisu[vue],forward,reverse,which=contenu_vue_,debut_vue=debut_vue,fin_vue=fin_vue,chr=chr,annot=annot,style=style,top=max_forward,bottom=max_reverse,marks=c(debut_vue,fin_vue),strandMarks=brin_transcrit,titres=description_data,repeated=TRUE,name_flags="flags",initialize_label_sizes=initialize_label_sizes)

	  cat("</div><hr/>",file=S,append=TRUE)

	  addTable(S,quantif[contenu_vue_,])
	  
	  closeHTML(S)
	}
      }
#       
      else {
	
        type_transcrits=unlist(strsplit(description_repeated$legend[i],split=","))
        n_type_transcrits=length(type_transcrits)
	##correction
	style_special=data.frame(type_transcrits,(1:n_type_transcrits)+1,rep(1,n_type_transcrits),rep("box",n_type_transcrits),rep(description_repeated$strand[i],n_type_transcrits),stringsAsFactors=FALSE)
	##
	names(style_special)=names(style)
	print(style_special)
	png(paste(repertoire,"/images/legend_",index_sequence[i],".png",sep=""),240,200)
	prev=par(no.readonly = TRUE)
	par(mar=c(0, 0, 2, 0)+0.1)
	plotLegend(style_special,100,100)
	par(prev)
	dev.off()

	jeton=TRUE
	big_gff = unique(rbind(annot,gff))
	gff_=big_gff[(big_gff$Chr==description_repeated$chr[i])&(big_gff$Start>=description_repeated$start[i])&(big_gff$Stop<=description_repeated$stop[i]),]
        for(j in type_transcrits) {
	  if(jeton) {
	    annot_special = gff_[gff_$Type==j,]
	    jeton=FALSE
	  }
	  else {
	    annot_special = rbind(annot_special,gff_[gff_$Type==j,])
	  }
	}
	annot_special=annot_special[sort(annot_special$Start,index.return=TRUE)$ix,]
	n_transcrit_special = dim(annot_special)[1]
	gene_=numeric(n_transcrit_special)
	for(j in 1:n_transcrit_special) {
	  tmp=unique(unlist(strsplit(c(annot_special[j,8],annot_special[j,9],annot_special[j,10]),",")))
	  tmp=tmp[tmp!="."]
	  gene_[j]=paste(tmp,collapse="<br/>")
	}
	SGD_special=paste("<a href=\"http://www.yeastgenome.org/cgi-bin/locus.fpl?locus=",annot_special$ID,"\"><img src=\"boutons/SGD_logo.gif\"/></a>",sep="")
	
        disp_annot_special=data.frame(annot_special[,c(6,2)],gene_,annot_special[,3:5],SGD_special,stringsAsFactors = FALSE)
	names(disp_annot_special)[3]="gene"
	names(disp_annot_special)[7]="SGD"

	jeton=TRUE
	for(j in 1:length(description_data)) {
	  density_ = numeric(n_transcrit_special)
	  if(data_repeated[[i]][[j]]$stranded) {
	    for(k in 1:n_transcrit_special) {
	      if(annot_special$Strand[k]=="+") {
			tmp=sum(data_repeated[[i]][[j]]$F[[1]][[index_norm]][(data_repeated[[i]][[j]]$F[[1]][[1]]<=(annot_special$Stop[k]-debut_vue+1))&(data_repeated[[i]][[j]]$F[[1]][[2]]>=(annot_special$Start[k]-debut_vue+1))])
	      }
	      else {
			if(annot_special$Strand[k]=="-") {
				tmp=sum(data_repeated[[i]][[j]]$R[[1]][[index_norm]][(data_repeated[[i]][[j]]$R[[1]][[1]]<=(annot_special$Stop[k]-debut_vue+1))&(data_repeated[[i]][[j]]$R[[1]][[2]]>=(annot_special$Start[k]-debut_vue+1))])
			}
			else {
				tmp=sum(data_repeated[[i]][[j]]$F[[1]][[index_norm]][(data_repeated[[i]][[j]]$F[[1]][[1]]<=(annot_special$Stop[k]-debut_vue+1))&(data_repeated[[i]][[j]]$F[[1]][[2]]>=(annot_special$Start[k]-debut_vue+1))])
				tmp=tmp+sum(data_repeated[[i]][[j]]$R[[1]][[index_norm]][(data_repeated[[i]][[j]]$R[[1]][[1]]<=(annot_special$Stop[k]-debut_vue+1))&(data_repeated[[i]][[j]]$R[[1]][[2]]>=(annot_special$Start[k]-debut_vue+1))])
			}
	      }
	      density_[k] = tmp/(annot_special$Stop[i]-annot_special$Start[i]+1)
	    }
	  }else {
	    for(k in 1:n_transcrit_special) {
	      tmp=sum(data_repeated[[i]][[j]]$F[[1]][[index_norm]][(data_repeated[[i]][[j]]$F[[1]][[1]]<=(annot_special$Stop[k]-debut_vue+1))&(data_repeated[[i]][[j]]$F[[1]][[2]]>=(annot_special$Start[k]-debut_vue+1))])
	      tmp=tmp+sum(data_repeated[[i]][[j]]$R[[1]][[index_norm]][(data_repeated[[i]][[j]]$R[[1]][[1]]<=(annot_special$Stop[k]-debut_vue+1))&(data_repeated[[i]][[j]]$R[[1]][[2]]>=(annot_special$Start[k]-debut_vue+1))])
	      density_[k] = tmp/(annot_special$Stop[i]-annot_special$Start[i]+1)
	    }
	  }
	  if(jeton) {
	    density_special=data.frame(density_,stringsAsFactors=FALSE)
	    jeton=FALSE
	  }else {
	    density_special=data.frame(density_special,density_,stringsAsFactors=FALSE)
	  }
	}
	names(density_special) = paste(description_data,"densities")
	N=dim(density_special)[2]
	for(j in 1:N) {
	  density_special[,j]=format(density_special[,j],scientific=TRUE,digits=3)
	  density_special[,j]=gsub("+","",density_special[,j],fixed=TRUE)
	  density_special[,j]=gsub("e0","e",density_special[,j],fixed=TRUE)
	  density_special[,j]=gsub("e-0","e-",density_special[,j],fixed=TRUE)
	  density_special[,j]=paste(gsub("e"," × 10<sup>",density_special[,j]),"</sup>",sep="")
	  density_special[,j]=gsub(" × 10<sup>0</sup>","",density_special[,j],fixed=TRUE)
	}

	for(vue in 1:n_vue) {
	  cat(paste("<?php \nheader('Location: ",name,"_",vue,".htm');\n ?>",sep=""),file=paste(repertoire,"/images/",toupper(name),"_",vue,".php",sep=""))
	  S=paste(repertoire,"/images/",name,"_",vue,".htm",sep="")
	  newHTML(S,paste(name,"  ",chr,":",debut_vue,"..",fin_vue," - ",description_vue[vue],sep=""))

	  addNavigatorHautRepeated(S,name,vue,description_vue,decal,l=paste("legend_",index_sequence[i],".png",sep=""))

	  contenu_vue_=contenu_vue[[vue]]
	  n_element_vue=length(contenu_vue_)
	  max_forward = numeric(n_element_vue)
	  max_reverse = numeric(n_element_vue)
	  for(element in 1:n_element_vue) {
	    i_data = contenu_vue_[element]
	    max_forward[element]=max(forward[[i_data]])
	    max_reverse[element]=max(reverse[[i_data]])
	  }
	  groupe_maximum_=groupe_maximum[[vue]]
	  for(i_max in unique(groupe_maximum_)) {
	    max_forward[groupe_maximum_==i_max]=max(max_forward[groupe_maximum_==i_max])
	    max_reverse[groupe_maximum_==i_max]=max(max_reverse[groupe_maximum_==i_max])
	  }
	  nom_fichier_image=paste(name,"_",vue,".png",sep="")
	  cat(paste("<p>\n<img src=\"",nom_fichier_image,"\" />\n</p>\n",sep=""),file=S,append=TRUE)
	  chemin_fichier_image=paste(repertoire,"/images/",nom_fichier_image,sep="")
	  
	  plotVisu(chemin_fichier_image,typeVisu=typeVisu[vue],forward,reverse,which=contenu_vue_,debut_vue=debut_vue,fin_vue=fin_vue,chr=chr,annot=annot_special,style=style_special,top=max_forward,bottom=max_reverse,strandMarks=brin_transcrit,titres=description_data,repeated=TRUE,name_flags="flags",display_name=display_name,initialize_label_sizes=initialize_label_sizes)

	  cat("</div><hr/>",file=S,append=TRUE)

	  annot_ = cbind(disp_annot_special,density_special[,contenu_vue_])
	  colnames(annot_) <- c(colnames(disp_annot_for_count), colnames(densitiesGenesNormal)[contenu_vue_])
	  addTable(S,annot_)
	  
	  closeHTML(S)
	}
      
      return(0)
    }
    }
    else {
      return(1)
    }
  }
  ctl_repeat=lapply(1:length(index_sequence), repeatProcess)
#     ctl_repeat=mclapply(1:length(index_sequence), repeatProcess, mc.cores=ncore)

}

cat("Done ! (check warnings)\n")
