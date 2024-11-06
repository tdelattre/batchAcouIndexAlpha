# ______________________________________________________________________________________ 
# ______________________________________________________________________________________ 
# ______________________________________________________________________________________ 
# |                                                                                    |
# |          SCRIPT WRITTEN BY THOMAS DELATTRE thomas.delattre@inrae.fr                | 
# |    99.9% Based on the work of Amandine Gasc! (see below and credit accordingly)    |
# |                             source material here :                                 |
# |   https://github.com/agasc/Soundscape-analysis-with-R/blob/master/AcouIndexAlpha.r |
# |                              ----------------                                      | 
# |                              LICENCE CC-BY-SA                                      | 
# |                              ----------------                                      |
# | This license lets others remix, adapt, and build upon your work even for           |
# | commercial purposes, as long as they credit you and license their new creations    |
# | under the identical terms.                                                         |
# |                                                                                    |
# | The proposed code has a purely academic purpose, is valid under the conditions     |
# | of use of the scientific project for which it was funded and at the date of        |
# | acceptance of the article presenting the code. As with any research work, the      |
# | code is not free of possible errors, approximations, sub-optimisations or          |
# | defects in monitoring dependencies between libraries of the program.               |
# |                                                                                    |
# ______________________________________________________________________________________ 
# |                                                                                    |
# | Cette licence permet à d'autres personnes de remixer, d'adapter et de              |
# | développer ce travail, même à des fins commerciales, à condition qu'elles          |
# | créditent l'auteur et accordent une licence pour leurs nouvelles créations aux     |
# | mêmes conditions.                                                                  |
# |                                                                                    |
# | Le code proposé a une visée purement académique, est valable dans les conditions   |
# | d'utilisation du projet scientifique pour lequel il a été financé et à la date de  |
# | d'acceptation de l'article de présentation du code.                                |
# | Comme tout travail de recherche, le code n'est pas exempt d'éventuelles erreurs,   |
# | approximations, sous-optimisations ou défauts de suivi des dépendances entre       |
# | sous-éléments du programme.                                                        |
# ______________________________________________________________________________________ 
# Objectif du script : appliquer les indices acoustiques sur les fichiers wav bruts d'une minutes collectés par des
# capteurs wildlife acoustics, sur la base d'une structure de fichiers
# projet / session / parcelle / Data2 / fichiers sons
# Le script consiste essentiellement en une sérialisation de la fonction AcouIndexAlphaCustom créée par Amandine Gasc,
# qui calcule les indices acoustiques les plus communs (ACI, NDSI, BI, ADI, NP)
# cette fonction est dans le fichier AcouIndexAlphaCustomSilent.r, merci de télécharger les deux fichiers
#-------------------------------------------------------------------------------------------------------------------
#changelog : 
## 28.05.24
# ajout de texte de convénience pour aller chercher les fichiers, depuis serveur ou machine perso
# ajout d'une barre de progression dans la fonction d'analyse acoustique, dans ldply (barre de progression interne à la parcelle)
# modification de AcouIndexAlphaCustom pour en faire une version "Silent", plus rapide car sans affichage des messages
# (cf. AcouIndexAlphaCustomSilent.R)
# 29.05 (v2) ajout du multitasking !! 
# 31.05 nettoyage et commentaires

#-------------------------------------------------------------------------------------------------------------------
library(stringr)
library(plyr)
library(parallel)
library(dplyr)
#-------------------------------------------------------------------------------------------------------------------
#Script source indices acoustiques (adapted from Amandine Gasc)
#---
source("/home/birdnet/AcouIndexAlphaCustomSilent.r")
#-------------------------------------------------------------------------------------------------------------------
#### Source des données acoustiques
#### DEPUIS SERVEUR : 
#rootdir="/home/partage/pshnas2/"
rootdir="/cbc/birdnet/pshnas3/"

#### DEPUIS POSTE PERSO : 
#rootdir="/home/tdelattre/pshnas2/"

#complément de chemin d'accès
rootdir=paste(rootdir,"DATA_2024/Automne_2024/",sep="")

#-------------------------------------------------------------------------------------------------------------------
#retirer X minutes (1 minute = 1 fichier wav) au début et fin pour éviter les bruits de manipulation et les fichiers foireux
cropBegin=15
cropEnd=15
#sessions (semaines) sur lesquelles appliquer le script
sessions=c("S1","S2")
#sessions=c("aig")
  #------------------------------------------ début de boucle
  #alldata=NULL #(ré)initialisation de l'objet pour éviter les bêtises
  
  #create an empty data frame with 0 rows and 10 columns
  alldata <- data.frame(matrix(ncol = 10, nrow = 0))
  #provide column names
  colnames(alldata) <- c("","Bioac_left","ACI_left","NDSI_left","ADI_left","npic_left","session","parcelle","mydir","fileName")
  alldata[1,]=rep(NA,length(alldata[1,]))
  for (j in 1:length(sessions)) {
    print("----------------------------------------------------------") #affichage à l'écran d'infos de base sur les données
    session_name=sessions[j]
    print(paste("session",
                session_name,
                "(",
                j,
                "/",
                length(sessions),
                ")"
                ))
    mydir = paste(rootdir,session_name,sep="")
    print(mydir)
    parcelles=list.dirs(path = mydir, full.names = FALSE, recursive = FALSE)
    print (paste("parcelles dans la session ", j, ":"))
    print(parcelles)
    
   for (i in 1:length(parcelles)) {
      #for (i in 1:12) {
      parcel_name=parcelles[i]
      print(paste("analyse de la parcelle",
                  parcel_name,
                  "(",
                  i,
                  "/",
                  length(parcelles),
                  ")")) #affiche à quelle parcelle nous en sommes
      #print("debug1")
      mydir = paste(rootdir,session_name,"/",parcel_name,"/Data2/",sep="") #à adapter en fonction du sous-dossier d'import
      myfiles = list.files(path=mydir, pattern="*.wav", full.names=FALSE) #valable tant que Wildlife acoustics respecte son système d'arborescence
      #print("debug2")
      # #lit tous les fichiers de la liste, y ajoute une colonne nom de fichier, les colle entre eux
      # dat_txt = ldply(myfiles, readCustom, sep = ",", fill=TRUE, header = TRUE,quote="")
      # 
      # #dat_txt$Begin.Time..s.=as.integer(dat_txt$Begin.Time..s.)#pour éviter les problèmes de format sur certains fichiers
      # 
      # #ajout de colonnes avec des informations fréquemment utiles (facultatif)
      # if (length(row.names(dat_txt)) > 0) { #évite les blocages en cas de fichiers vides
      #   dat_txt$session=session_name #ajoute l'information session et parcelle, disponibles uniquement à ce moment
      #   dat_txt$parcelle=parcel_name
      #   dat_txt$mydir=mydir
      # }
      
      #lit les fichiers wav de la liste (sauf X premiers et derniers) en calcule les indices, combine les résultats
      #dat_txt = ldply(myfiles[cropBegin:(length(myfiles)-cropEnd)], AcouIndexAlpha, mydir=mydir, min_freq = 2000, max_freq = 12000,.progress="text")
      
      #version calcul parrallèle du ldply
      #--------------- début du parallel processing
      cl=makeCluster(detectCores()-34) #création d'un cluster de processeurs prenant tous les coeurs dispos, -4
      
      #print("debug2.1")
      #export des fonctions "maison" dans le cluster de processeurs
      clusterExport(cl, c("AcouIndexAlpha","bioacSilent",
                          "acoustic_complexity_Silent",
                          "ndsi_Silent","acoustic_diversity_Silent",
                          "fpeaksFlat"), 
                    envir=environment()) 
      #print("debug2.2")
      print(paste("processing",
                  length(myfiles[cropBegin:(length(myfiles)-cropEnd)]),
                  "files with a socket cluster containing", 
                  length(cl), 
                  "cores. Processing started at",
                  Sys.time()))
      #calcul en parallèle des indices acoustiques en utilisant AcouIndexAlpha
      dat_txt = parLapply(cl,
                          myfiles[cropBegin:(length(myfiles)-cropEnd)], 
                          AcouIndexAlpha,
                          mydir=mydir,
                          min_freq = 2000,
                          max_freq = 12000) 
      #print("debug2.3")
      stopCluster(cl) 
      #--------------- fin du parallel processing
      
      #parLapply crée une liste de dataframes : on les aplatit en un seul df
      dat_txt=bind_rows(dat_txt)
      
      #print("debug3")
      # #ajout de colonnes avec des informations fréquemment utiles (facultatif)
      #if (length(row.names(dat_txt)) > 0) { #évite les blocages en cas de fichiers vides
      if (length(dat_txt) > 0) { #évite les blocages en cas de fichiers vides
        dat_txt$session=session_name #ajoute l'information session et parcelle, disponibles uniquement à ce moment
        dat_txt$parcelle=parcel_name
        dat_txt$mydir=mydir
        dat_txt$fileName=myfiles[cropBegin:(length(myfiles)-cropEnd)]
      }
      #print("debug4")
      alldata=dplyr::bind_rows(alldata,dat_txt) #colle les sorties des différentes parcelles et sessions en un objet
      #print("debug5")
    }
  }
  #------------------------------------------ Fin de boucle
  
  #on écrit alldata (compilation) sauf ligne 1 (NAs)
  write.csv2(alldata[2:length(alldata$Bioac_left),],paste(rootdir,"5_acoustic_index_",str_flatten(sessions),"_2emeTentative.txt",sep=""))
  
  #exploration rapide des données
  #GGally::ggpairs(alldata)

