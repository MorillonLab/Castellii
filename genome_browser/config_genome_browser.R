#  ********************
#  *  load functions  * 
#  ********************

source("/home/ugo/SCRIPTS/BAMTOHTML/bamToHTML_functions_v2.6.1.R")

#  ***************
#  *  variables  *
#  ***************

# string, title html index page

titre_page=""

# string, folder name

texte_description="."

# string, folder of files used in main, style file, design.css, ...
common_files_dir : répertoire des fichiers style, design.css, etc...

# integer, size visualisation window 

taille_fenetre = 15000

# integer, space between windows

decal = 5000

# Data frame, Type (transcripts type, as in annotation file), col (color), label (ID shown)

style = read.table(paste0(common_files_dir, "style.tab"), sep="\t",as.is=TRUE,header=TRUE)

# Data frame, obtained through readGff function

annot = readGff("annot.gff")

# vector, transcripts type not linked to SGD

nonSGD = c("")	

# integer vector, from 1 to number of bamHandler in data

index_data = c()

# sample names, length should be equal to index_data

description_data = c("")

# string vector, bam files (sorted, indexed)

file_names = c("")

# integer, number of core used in computation

ncore = 1

# bamHandler list, obtained from readBam function
# input : file : bam files
#         stranded = T or F (strand specific or not)
#         libraryType = "standard" or "inverse" (if stranded is True)
#         
# output : list of bamHandler, 1 per bam file 

datarough = readBam(file = file_names, stranded = T/F , libraryType=c("standard","inverse"))

# numeric vector, gives normalization coefficient

poids = c()

# normalized bamHandler

data=normalizeBams(datarough, poids, relative = F)

# to pool samples from the same condition, one can use the function meanBam
# for exemple, in data at index k will be stored the mean signal of the content of data at indices i and j
# data[[k]] = meanBam(data[[i]], data[[j]])

# list of vector, each containing indexes of data element

contenu_vue = list()

# string vector, description of contenue_vue 

description_vue = c("")

# string vector, type of visualisation for each group ("heatmap", "lines" or "classic")

typeVisu = c("")  

# boolean, should signal be log transformed (log2)

isVisuLog = T or F

# launch genome browser

source("~/SCRIPTS/BAMTOHTML/bamToHTML_main_v2.6_corrected.R") 
