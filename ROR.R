##################################
#
#   analyse donnees vacci CD 2A -- ror
#
# V1  août 2017
##################################




rm(list=ls())
setwd("I:\\DIRECTION ACTION TERRITORIALE DE SANTE\\VSSE\\CIRE\\Dossiers MI\\Vaccination\\etude_CG2A\\analyse")


load("donnees/df_pour_analyse.RData")

ror <- c("O","ROR","ROU","RRUB","RUB")


# interet a separer les 3 valences 
tempo <- donnees[donnees$vaccin_code %in% ror,] 
table(factor(tempo$vaccin_code)) # affichage de la repartition
rm(tempo) # on fait de la place

# => aucun interet à separer, car pourcentage très faible, on restreint aux vacci avec les 3 valences
ror <- c("ROR")

####################################################################
# analyse donnees nationales
###################################################################



# fin analyse nationale ror     
#################################

#####################################################################
# creation de la fonction calcul CV 24 mois pour tempo_ror
f_ror24 <- function (df_ror,an){
      ############################################################
      # preparation
      # df_ror <- donnees ; an <- 1997
      # restriction du df d'entree a l'annee voulue et au vaccin
      ddn_min <- as.Date(paste(an-2,"-01-01",sep="")) # on veut les enfants ayants 24 mois 
      ddn_max <- as.Date(paste(an-2,"-12-31",sep=""))
      tempo_ror <- df_ror[df_ror$datenaiss > ddn_min &
                          df_ror$datenaiss < ddn_max,] # temp_ror ne contient que les enfants de 24 mois l'annee choisie
      
      tempo_population <- length(unique(tempo_ror$nodossier)) # calcul population denominateur avant restriction sur vaccins et age a la vacci
      
      # on garde les vacci effectuees avant 24 mois (2 ans = 730 jours) et vacci ror
      tempo_ror <- tempo_ror[tempo_ror$dateopv < (tempo_ror$datenaiss + 730)  &
                                   tempo_ror$vaccin_code %in% ror,]
      

      # on classe suivant le n° de dossier et l'anciennete de l'opv
      tempo_ror <- tempo_ror[order(tempo_ror$nodossier,tempo_ror$dateopv),]
      
      #############################################################
      # debut des calculs des nombres de vaccins
      
      # creation d'un df avec une seule repetition pour chaque dossier => nb de vacci 1 dose
      tempo_ror_1dose <- tempo_ror[!duplicated(tempo_ror$nodossier),]
      # pour 2 doses, on prend les autres vacci donc, sans "!", puis on réapplique
      # utilisation d'un df "pivot" tempo_ror_dose
      tempo_ror_dose <- tempo_ror[duplicated(tempo_ror$nodossier),] # toutes les vacci sauf la premiere dose
      tempo_ror_2dose <- tempo_ror_dose[!duplicated(tempo_ror_dose$nodossier),]
      
      
      sortie1 <- list(
            round(dim(tempo_ror_1dose)[1] / tempo_population * 100,2),
            binom.test(dim(tempo_ror_1dose)[1],tempo_population)$conf.int,
            summary(as.numeric(tempo_ror_1dose$age_vacci_mois))
            )
      sortie2 <- list(
            round(dim(tempo_ror_2dose)[1] / tempo_population * 100,2),
            binom.test(dim(tempo_ror_2dose)[1],tempo_population)$conf.int,
            summary(as.numeric(tempo_ror_2dose$age_vacci_mois))
            )
      
      sortie <- list(sortie1,sortie2) # liste de listes ...
      return(sortie)
      
}
########################

#####################################################################
# creation de la fonction calcul CV 4 ans pour tempo_ror
f_ror4 <- function (df_ror,an){
      ############################################################
      # preparation
      
      # restriction du df d'entree a l'annee voulue et au vaccin
      ddn_min <- as.Date(paste(an-4,"-01-01",sep="")) # on veut les enfants ayants 4 ans l'annee choisie
      ddn_max <- as.Date(paste(an-4,"-12-31",sep=""))
      tempo_ror <- df_ror[df_ror$datenaiss > ddn_min &
                          df_ror$datenaiss < ddn_max,]
      
      tempo_population <- length(unique(tempo_ror$nodossier)) # calcul population denominateur avant restriction sur vaccins
      
      # on garde les vacci effectuees avant 4 ans (= 1460 jours) et sur ror
      tempo_ror <- tempo_ror[tempo_ror$dateopv < (tempo_ror$datenaiss + 1460) &
                                   tempo_ror$vaccin_code %in% ror,]
      
      
      # on classe suivant le n° de dossier et l'anciennete de l'opv
      tempo_ror <- tempo_ror[order(tempo_ror$nodossier,tempo_ror$dateopv),]
      
      #############################################################
      # debut des calculs des nombres de vaccines
      
      # creation d'un df avec une seule repetition pour chaque dossier => nb de vacci 1 dose
      tempo_ror_1dose <- tempo_ror[!duplicated(tempo_ror$nodossier),]
      # pour 2 doses, on prend les autres vacci donc, sans "!", puis on réapplique
      # utilisation d'un df "pivot" tempo_ror_dose
      tempo_ror_dose <- tempo_ror[duplicated(tempo_ror$nodossier),] # toutes les vacci sauf la premiere dose
      tempo_ror_2dose <- tempo_ror_dose[!duplicated(tempo_ror_dose$nodossier),]
      
      
      sortie1 <- list(
            round(dim(tempo_ror_1dose)[1] / tempo_population * 100,2),
            binom.test(dim(tempo_ror_1dose)[1],tempo_population)$conf.int,
            summary(as.numeric(tempo_ror_1dose$age_vacci_mois))
            )
      sortie2 <- list(
            round(dim(tempo_ror_2dose)[1] / tempo_population * 100,2),
            binom.test(dim(tempo_ror_2dose)[1],tempo_population)$conf.int,
            summary(as.numeric(tempo_ror_2dose$age_vacci_mois))
            )
      
      sortie <- list(sortie1,sortie2) # liste de listes ...
      return(sortie)
}
########################
###################################################
# boucle de sortie des resultats

sink("sorties/ror/ror_24.txt") # sorties vers le fichier specifie

print("###################################### 24 mois #########################################")

#remplissage des tableaux de sortie
for (d in 1:2){ # d pour dose
      print(paste("---------------------------",d," dose(s) ----------------------------"))
      
      sortie <- matrix(NA,7,19)
      colnames(sortie) <- 1995:2013
      rownames(sortie) <- c("CV en %","IC inf (%)","IC sup (%)","1er quartile","médiane","moyenne","3ième quartile")
   
      for (i in 1995:2013) {
            tempo <- f_ror24(donnees,i)[[d]] # on ne recupere que la "d"ieme liste dans la sortie
            sortie[1,i-1994] <- tempo[[1]] # CV
            sortie[2,i-1994] <- tempo[[2]][1]*100 ; sortie[3,i-1994] <- tempo[[2]][2]*100 # IC
            for (j in 4:7)
                  sortie[j,i-1994] <- tempo[[3]][j-2] # repartition des ages en mois
      }
      print(round(sortie,2))
      rm(tempo,sortie) # on efface les fichiers temporaires pour faire de la place
}

sink()

sink("sorties/ror/ror_4.txt") # sorties vers le fichier specifie

print("###################################### 4 ans #########################################")

#remplissage des tableaux de sortie
for (d in 1:2){ # d pour dose
      print(paste("---------------------------",d," dose(s) ----------------------------"))
      
      sortie <- matrix(NA,7,19)
      colnames(sortie) <- 1997:2015
      rownames(sortie) <- c("CV en %","IC inf (%)","IC sup (%)","1er quartile","médiane","moyenne","3ième quartile")
   
      for (i in 1997:2015) {
            tempo <- f_ror4(donnees,i)[[d]] # on ne recupere que la "d"ieme liste dans la sortie
            sortie[1,i-1996] <- tempo[[1]] # CV
            sortie[2,i-1996] <- tempo[[2]][1]*100 ; sortie[3,i-1996] <- tempo[[2]][2]*100 # IC
            for (j in 4:7)
                  sortie[j,i-1996] <- tempo[[3]][j-2] # repartition des ages en mois
      }
      print(round(sortie,2))
      rm(tempo,sortie) # on efface les fichiers temporaires pour faire de la place
}

sink()

sink("sorties/ror/ror_canton.txt") # sorties vers le fichier specifie


cat("
      ******************************
      ****** ror *******************
      ******************************
      ****** PAR CANTONS ***********
      ******************************")

#####################################################################
# fonction calcul CV 24 mois pour tempo_ror pour les cantons
f_ror24_canton <- function (df_ror){
      ############################################################
      # preparation
      # df_ror <- donnees[donnees$,] ; an <- 1997
      # restriction du df d'entree a l'annee voulue et au vaccin
      ddn_min <- as.Date("2007-01-01") # on veut les enfants ayants 24 mois 
      ddn_max <- as.Date("2011-12-31")
      tempo_ror <- df_ror[df_ror$datenaiss > ddn_min &
                          df_ror$datenaiss < ddn_max,] # temp_ror ne contient que les enfants de 24 mois l'annee choisie
      
      # calcul population denominateur avant restriction sur vaccins et age
      tempo_population <- length(unique(tempo_ror$nodossier))
      
      # on garde les vacci effectuees avant 24 mois (2 ans = 730 jours) et en lien ror
      tempo_ror <- tempo_ror[tempo_ror$dateopv < (tempo_ror$datenaiss + 730) &
                                   tempo_ror$vaccin_code %in% ror,]

      # on classe suivant le n° de dossier et l'anciennete de l'opv
      tempo_ror <- tempo_ror[order(tempo_ror$nodossier,tempo_ror$dateopv),]
      
      #############################################################
      # debut des calculs des nombres de vaccins
      
      # creation d'un df avec une seule repetition pour chaque dossier => nb de vacci 1 dose
      tempo_ror_1dose <- tempo_ror[!duplicated(tempo_ror$nodossier),]
      # pour 2 doses, on prend les autres vacci donc, sans "!", puis on réapplique
      # utilisation d'un df "pivot" tempo_ror_dose
      tempo_ror_dose <- tempo_ror[duplicated(tempo_ror$nodossier),] # toutes les vacci sauf la premiere dose
      tempo_ror_2dose <- tempo_ror_dose[!duplicated(tempo_ror_dose$nodossier),]
      
      sortie <- cat(
            "vaccination ror dans le canton à 24 mois entre 2009 et 2013 en % : ",
            round(length(unique(tempo_ror_1dose$nodossier))/tempo_population*100,1),"\n",
            "intervalle de confiance : ",
            round(binom.test(dim(tempo_ror_1dose)[1],tempo_population)$conf.int*100,1),"\n",
            "2ieme vaccination ror dans le canton à 24 mois entre 2009 et 2013 en % : ",
            round(length(unique(tempo_ror_2dose$nodossier))/tempo_population*100,1),"\n",
            "intervalle de confiance : ",
            round(binom.test(dim(tempo_ror_2dose)[1],tempo_population)$conf.int*100,1)
      )
      return(sortie)
      
}
########################


#####################################################################
# fonction calcul CV 4 ans pour tempo_ror pour les cantons
f_ror4_canton <- function (df_ror){
      ############################################################
      # preparation
      # df_ror <- donnees[donnees$,] ; an <- 1997
      # restriction du df d'entree a l'annee voulue et au vaccin
      ddn_min <- as.Date("2007-01-01") # on veut les enfants ayants 24 mois 
      ddn_max <- as.Date("2011-12-31")
      tempo_ror <- df_ror[df_ror$datenaiss > ddn_min &
                          df_ror$datenaiss < ddn_max,] # temp_ror ne contient que les enfants de 24 mois l'annee choisie
      
      # calcul population denominateur avant restriction sur vaccins et age
      tempo_population <- length(unique(tempo_ror$nodossier))
      
      # on garde les vacci effectuees avant 4 ans (4 ans = 1460 jours) et en lien ror
      tempo_ror <- tempo_ror[tempo_ror$dateopv < (tempo_ror$datenaiss + 1460) &
                                   tempo_ror$vaccin_code %in% ror,]

      # on classe suivant le n° de dossier et l'anciennete de l'opv
      tempo_ror <- tempo_ror[order(tempo_ror$nodossier,tempo_ror$dateopv),]
      
      #############################################################
      # debut des calculs des nombres de vaccins
      
      # creation d'un df avec une seule repetition pour chaque dossier => nb de vacci 1 dose
      tempo_ror_1dose <- tempo_ror[!duplicated(tempo_ror$nodossier),]
      # pour 2 doses, on prend les autres vacci donc, sans "!", puis on réapplique
      # utilisation d'un df "pivot" tempo_ror_dose
      tempo_ror_dose <- tempo_ror[duplicated(tempo_ror$nodossier),] # toutes les vacci sauf la premiere dose
      tempo_ror_2dose <- tempo_ror_dose[!duplicated(tempo_ror_dose$nodossier),]
      
      sortie <- cat(
            "vaccination ror dans le canton à 4 ans entre 2011 et 2015 en % : ",
            round(length(unique(tempo_ror_1dose$nodossier))/tempo_population*100,1),"\n",
            "intervalle de confiance : ",
            round(binom.test(dim(tempo_ror_1dose)[1],tempo_population)$conf.int*100,1),"\n",
            "2ième vaccination ror dans le canton à 4 ans entre 2011 et 2015 en % : ",
            round(length(unique(tempo_ror_2dose$nodossier))/tempo_population*100,1),"\n",
            "intervalle de confiance : ",
            round(binom.test(dim(tempo_ror_2dose)[1],tempo_population)$conf.int*100,1)
      )
      return(sortie)
      
}

######################################################################


for (c in levels(donnees$canton)){
      cat("
          ************************",c,"*********************************")
      
      tempo_canton <- donnees[donnees$canton == c,]
      cat("\n")
      f_ror24_canton(tempo_canton)
      cat("\n")
      f_ror4_canton(tempo_canton)
}


# ###### analyse des donnees


moyenne <- read.csv2("sorties/ror/ror24_3_canton_pour moyenne.csv")

round(addmargins(moyenne,2,mean)[,6],1)

result <- matrix(round(as.numeric(addmargins(moyenne,2,mean)[,6]),1)[c(T,F,F)],2,6)

rownames(result) <- c("vacci","2ième vacci")
colnames(result) <- levels(donnees$canton)
write.csv2(result,"sorties/ror/result_cantons.csv")

sink()

cat("
      **********************************
      **********************************
      ****** ANALYSE AGE PRIMOVACCI ****
      **********************************
      **********************************")


      vacci <- donnees[donnees$vaccin_code %in% ror,]

      # on classe suivant le n° de dossier et l'anciennete de l'opv
      vacci <- vacci[order(vacci$nodossier,vacci$dateopv),]
      
      #############################################################
      
      # contient 1 vacci
      Prem_vacci <- vacci[!duplicated(vacci$nodossier),] # seulement la premiere vacci
      
      # rappelvacci <- primovacci     # en attente, pour analyse age au rappel
      
      # analyse en tant que telle
      summary(Prem_vacci$age_vacci_mois)
      # analyse de ce df
      png("sorties/ror/age-Prem_vacci%02d.png", width=25 ,height=15, units="cm",res = 400)
            barplot(table(Prem_vacci$age_vacci_mois[Prem_vacci$age_vacci_mois < 48]),
                 main = "Répartition de l'âge en mois à la première vaccination ROR,
                  pour les enfants nés entre le 1er janvier 1993
                  et le 31 décembre 2011, Corse-du-Sud",
                 xlab = "âge en mois - coupure à 48 mois", ylab = "effectif") 
      dev.off()

      

cat("
      **********************************
      **********************************
      ******** ANALYSE AGE 2ieme vaccination ******
      **********************************
      **********************************")


      #############################################################
      # on reapplique, donc, 2 vacci
      Deuz_vacci <- vacci[duplicated(vacci$nodossier),]
      # mais pas plus, pour avoir l'age au dernier acte
      Deuz_vacci <- Deuz_vacci[!duplicated(Deuz_vacci$nodossier),]      

      # analyse en tant que telle
      summary(Deuz_vacci$age_vacci_mois)
      # analyse de ce df
      png("sorties/ror/age-deuxieme_vacci%02d.png", width=25 ,height=15, units="cm",res = 400)
            barplot(table(Deuz_vacci$age_vacci_mois[Deuz_vacci$age_vacci_mois < 48]),
                  main = "Répartition de l'âge en mois lors de la deuxième vaccination ROR,
                  pour les enfants nés entre le 1er janvier 1993
                  et le 31 décembre 2011, Corse-du-Sud",
                  xlab = "âge en mois - coupure à 48 mois", ylab = "effectif") 
      dev.off()





# fin analyse ror
####################################################################################################


#################################
#
# THAT'S ALL FOLKS !!!!!!!!!!!
#
#################################
