##################################
#
#   analyse donnees vacci CD 2A -- hib
#
# V1  août 2017
##################################

rm(list = ls())

setwd("I:\\DIRECTION ACTION TERRITORIALE DE SANTE\\VSSE\\CIRE\\Dossiers MI\\Vaccination\\etude_CG2A\\analyse")
load("donnees/df_pour_analyse.RData")

hib <- c("DTCPHIB","DTCPHIBHB","DTP CA+HIB","HIB")

tempo <- donnees[donnees$vaccin_code %in% hib,] 
table(factor(tempo$vaccin_code)) # affichage de la repartition
rm(tempo) # on fait de la place


####################################################################
# analyse donnees nationales
###################################################################

donnees_nat <- read.csv2("donnees/donnees_nationales_CS_24_hib.csv",header=F)

colnames(donnees_nat) <- 1998:2015
rownames(donnees_nat) <- c("3 doses")


# analyse evolution hib

      print(summary(lm(as.numeric(donnees_nat) ~ as.numeric(1998:2015))))

# # graphe evolution
# couleurs <- c("red")
# 
# png("sorties/hib/evolution_nat_CS24.png",width=25 ,height=15, units="cm",res = 400)
# matplot(t(donnees_nat),type="b",pch=1,xlab = "années", ylim = c(0,100), ylab = "taux de couverture vaccinal en %",
#          main = "Couverture vaccinale hib à l'âge de 24 mois, 
#          analyse des CS 24 - France, 1998-2015 - résultats en %", xaxt="n",
#         col = couleurs)
# axis(1,at=1:18,labels=1998:2015)
# abline(h = 95, col = "darkorange",lty = 2)
# dev.off()
 

# fin analyse nationale hib     
#################################

#####################################################################
# creation de la fonction calcul CV 24 mois pour tempo_hib
f_hib24 <- function (df_hib,an){
      ############################################################
      # preparation
      # df_hib <- donnees ; an <- 1997
      # restriction du df d'entree a l'annee voulue et au vaccin
      ddn_min <- as.Date(paste(an-2,"-01-01",sep="")) # on veut les enfants ayants 24 mois 
      ddn_max <- as.Date(paste(an-2,"-12-31",sep=""))
      tempo_hib <- df_hib[df_hib$datenaiss > ddn_min &
                          df_hib$datenaiss < ddn_max,] # temp_hib ne contient que les enfants de 24 mois l'annee choisie
      
      tempo_population <- length(unique(tempo_hib$nodossier)) # calcul population denominateur avant restriction sur vaccins et age a la vacci
      
      # on garde les vacci effectuees avant 24 mois (2 ans = 730 jours) et vacci hib
      tempo_hib <- tempo_hib[tempo_hib$dateopv < (tempo_hib$datenaiss + 730)  &
                                   tempo_hib$vaccin_code %in% hib,]
      

      # on classe suivant le n° de dossier et l'anciennete de l'opv
      tempo_hib <- tempo_hib[order(tempo_hib$nodossier,tempo_hib$dateopv),]
      
      #############################################################
      # debut des calculs des nombres de vaccins
      
      # creation d'un df avec une seule repetition pour chaque dossier => nb de vacci 1 dose
      tempo_hib_1dose <- tempo_hib[!duplicated(tempo_hib$nodossier),]
      # pour 2 doses, on prend les autres vacci donc, sans "!", puis on réapplique
      # utilisation d'un df "pivot" tempo_hib_dose
      tempo_hib_dose <- tempo_hib[duplicated(tempo_hib$nodossier),] # toutes les vacci sauf la premiere dose
      tempo_hib_2dose <- tempo_hib_dose[!duplicated(tempo_hib_dose$nodossier),]
      # idem 3ieme dose
      tempo_hib_dose <- tempo_hib_dose[duplicated(tempo_hib_dose$nodossier),]
      tempo_hib_3dose <- tempo_hib_dose[!duplicated(tempo_hib_dose$nodossier),]
      # idem 4ieme dose
      tempo_hib_dose <- tempo_hib_dose[duplicated(tempo_hib_dose$nodossier),]
      tempo_hib_4dose <- tempo_hib_dose[!duplicated(tempo_hib_dose$nodossier),]
      
      
      
      sortie1 <- list(
            round(dim(tempo_hib_1dose)[1] / tempo_population * 100,2),
            binom.test(dim(tempo_hib_1dose)[1],tempo_population)$conf.int,
            summary(as.numeric(tempo_hib_1dose$age_vacci_mois))
            )
      sortie2 <- list(
            round(dim(tempo_hib_2dose)[1] / tempo_population * 100,2),
            binom.test(dim(tempo_hib_2dose)[1],tempo_population)$conf.int,
            summary(as.numeric(tempo_hib_2dose$age_vacci_mois))
            )
      sortie3 <- list(
            round(dim(tempo_hib_3dose)[1] / tempo_population * 100,2),
            binom.test(dim(tempo_hib_3dose)[1],tempo_population)$conf.int,
            summary(as.numeric(tempo_hib_3dose$age_vacci_mois))
            )
      sortie4 <- list(
            round(dim(tempo_hib_4dose)[1] / tempo_population * 100,2),
            binom.test(dim(tempo_hib_4dose)[1],tempo_population)$conf.int,
            summary(as.numeric(tempo_hib_4dose$age_vacci_mois))
            )
      sortie <- list(sortie1,sortie2,sortie3,sortie4) # liste de listes ...
      return(sortie)
      
}
########################

#####################################################################
# creation de la fonction calcul CV 4 ans pour tempo_hib
f_hib4 <- function (df_hib,an){
      ############################################################
      # preparation
      
      # restriction du df d'entree a l'annee voulue et au vaccin
      ddn_min <- as.Date(paste(an-4,"-01-01",sep="")) # on veut les enfants ayants 4 ans l'annee choisie
      ddn_max <- as.Date(paste(an-4,"-12-31",sep=""))
      tempo_hib <- df_hib[df_hib$datenaiss > ddn_min &
                          df_hib$datenaiss < ddn_max,]
      
      tempo_population <- length(unique(tempo_hib$nodossier)) # calcul population denominateur avant restriction sur vaccins
      
      # on garde les vacci effectuees avant 4 ans (= 1460 jours) et sur hib
      tempo_hib <- tempo_hib[tempo_hib$dateopv < (tempo_hib$datenaiss + 1460) &
                                   tempo_hib$vaccin_code %in% hib,]
      
      
      # on classe suivant le n° de dossier et l'anciennete de l'opv
      tempo_hib <- tempo_hib[order(tempo_hib$nodossier,tempo_hib$dateopv),]
      
      #############################################################
      # debut des calculs des nombres de vaccines
      
      # creation d'un df avec une seule repetition pour chaque dossier => nb de vacci 1 dose
      tempo_hib_1dose <- tempo_hib[!duplicated(tempo_hib$nodossier),]
      # pour 2 doses, on prend les autres vacci donc, sans "!", puis on réapplique
      # utilisation d'un df "pivot" tempo_hib_dose
      tempo_hib_dose <- tempo_hib[duplicated(tempo_hib$nodossier),] # toutes les vacci sauf la premiere dose
      tempo_hib_2dose <- tempo_hib_dose[!duplicated(tempo_hib_dose$nodossier),]
      # idem 3ieme dose
      tempo_hib_dose <- tempo_hib_dose[duplicated(tempo_hib_dose$nodossier),]
      tempo_hib_3dose <- tempo_hib_dose[!duplicated(tempo_hib_dose$nodossier),]
      # idem 4ieme dose
      tempo_hib_dose <- tempo_hib_dose[duplicated(tempo_hib_dose$nodossier),]
      tempo_hib_4dose <- tempo_hib_dose[!duplicated(tempo_hib_dose$nodossier),]
      
      sortie1 <- list(
            round(dim(tempo_hib_1dose)[1] / tempo_population * 100,2),
            binom.test(dim(tempo_hib_1dose)[1],tempo_population)$conf.int,
            summary(as.numeric(tempo_hib_1dose$age_vacci_mois))
            )
      sortie2 <- list(
            round(dim(tempo_hib_2dose)[1] / tempo_population * 100,2),
            binom.test(dim(tempo_hib_2dose)[1],tempo_population)$conf.int,
            summary(as.numeric(tempo_hib_2dose$age_vacci_mois))
            )
      sortie3 <- list(
            round(dim(tempo_hib_3dose)[1] / tempo_population * 100,2),
            binom.test(dim(tempo_hib_3dose)[1],tempo_population)$conf.int,
            summary(as.numeric(tempo_hib_3dose$age_vacci_mois))
            )
      sortie4 <- list(
            round(dim(tempo_hib_4dose)[1] / tempo_population * 100,2),
            binom.test(dim(tempo_hib_4dose)[1],tempo_population)$conf.int,
            summary(as.numeric(tempo_hib_4dose$age_vacci_mois))
            )
      sortie <- list(sortie1,sortie2,sortie3,sortie4) # liste de listes ...
      return(sortie)
}
########################
###################################################
# boucle de sortie des resultats

sink("sorties/hib/hib_24.txt") # sorties vers le fichier specifie

print("###################################### 24 mois #########################################")

#remplissage des tableaux de sortie
for (d in 3:4){ # d pour dose
      print(paste("---------------------------",d," dose(s) ----------------------------"))
      
      sortie <- matrix(NA,7,19)
      colnames(sortie) <- 1995:2013
      rownames(sortie) <- c("CV en %","IC inf (%)","IC sup (%)","1er quartile","médiane","moyenne","3ième quartile")
   
      for (i in 1995:2013) {
            tempo <- f_hib24(donnees,i)[[d]] # on ne recupere que la "d"ieme liste dans la sortie
            sortie[1,i-1994] <- tempo[[1]] # CV
            sortie[2,i-1994] <- tempo[[2]][1]*100 ; sortie[3,i-1994] <- tempo[[2]][2]*100 # IC
            for (j in 4:7)
                  sortie[j,i-1994] <- tempo[[3]][j-2] # repartition des ages en mois
      }
      print(round(sortie,2))
      rm(tempo,sortie) # on efface les fichiers temporaires pour faire de la place
}

sink()

sink("sorties/hib/hib_4.txt") # sorties vers le fichier specifie
print("###################################### 4 ans #########################################")

#remplissage des tableaux de sortie
for (d in 3:4){ # d pour dose
      print(paste("---------------------------",d," dose(s) ----------------------------"))
      
      sortie <- matrix(NA,7,19)
      colnames(sortie) <- 1997:2015
      rownames(sortie) <- c("CV en %","IC inf (%)","IC sup (%)","1er quartile","médiane","moyenne","3ième quartile")
   
      for (i in 1997:2015) {
            tempo <- f_hib4(donnees,i)[[d]] # on ne recupere que la "d"ieme liste dans la sortie
            sortie[1,i-1996] <- tempo[[1]] # CV
            sortie[2,i-1996] <- tempo[[2]][1]*100 ; sortie[3,i-1996] <- tempo[[2]][2]*100 # IC
            for (j in 4:7)
                  sortie[j,i-1996] <- tempo[[3]][j-2] # repartition des ages en mois
      }
      print(round(sortie,2))
      rm(tempo,sortie) # on efface les fichiers temporaires pour faire de la place
}

sink()

sink("sorties/hib/hib_canton.txt") # sorties vers le fichier specifie
cat("
      ******************************
      ****** hib *******************
      ******************************
      ****** PAR CANTONS ***********
      ******************************")

#####################################################################
# fonction calcul CV 24 mois pour tempo_hib pour les cantons
f_hib24_canton <- function (df_hib){
      ############################################################
      # preparation
      # df_hib <- donnees[donnees$,] ; an <- 1997
      # restriction du df d'entree a l'annee voulue et au vaccin
      ddn_min <- as.Date("2007-01-01") # on veut les enfants ayants 24 mois 
      ddn_max <- as.Date("2011-12-31")
      tempo_hib <- df_hib[df_hib$datenaiss > ddn_min &
                          df_hib$datenaiss < ddn_max,] # temp_hib ne contient que les enfants de 24 mois l'annee choisie
      
      # calcul population denominateur avant restriction sur vaccins et age
      tempo_population <- length(unique(tempo_hib$nodossier))
      
      # on garde les vacci effectuees avant 24 mois (2 ans = 730 jours) et en lien hib
      tempo_hib <- tempo_hib[tempo_hib$dateopv < (tempo_hib$datenaiss + 730) &
                                   tempo_hib$vaccin_code %in% hib,]

      # on classe suivant le n° de dossier et l'anciennete de l'opv
      tempo_hib <- tempo_hib[order(tempo_hib$nodossier,tempo_hib$dateopv),]
      
      #############################################################
      # debut des calculs des nombres de vaccins
      
      # creation d'un df avec une seule repetition pour chaque dossier => nb de vacci 1 dose
      tempo_hib_1dose <- tempo_hib[!duplicated(tempo_hib$nodossier),]
      # pour 2 doses, on prend les autres vacci donc, sans "!", puis on réapplique
      # utilisation d'un df "pivot" tempo_hib_dose
      tempo_hib_dose <- tempo_hib[duplicated(tempo_hib$nodossier),] # toutes les vacci sauf la premiere dose
      tempo_hib_2dose <- tempo_hib_dose[!duplicated(tempo_hib_dose$nodossier),]
      # idem 3ieme dose
      tempo_hib_dose <- tempo_hib_dose[duplicated(tempo_hib_dose$nodossier),]
      tempo_hib_3dose <- tempo_hib_dose[!duplicated(tempo_hib_dose$nodossier),]
      # idem 4ieme dose
      tempo_hib_dose <- tempo_hib_dose[duplicated(tempo_hib_dose$nodossier),]
      tempo_hib_4dose <- tempo_hib_dose[!duplicated(tempo_hib_dose$nodossier),]
      
      sortie <- cat(
            "primovaccination hib dans le canton à 24 mois entre 2009 et 2013 en % : ",
            round(length(unique(tempo_hib_3dose$nodossier))/tempo_population*100,1),"\n",
            "intervalle de confiance : ",
            round(binom.test(dim(tempo_hib_3dose)[1],tempo_population)$conf.int*100,1),"\n",
            "rappel hib dans le canton à 24 mois entre 2009 et 2013 en % : ",
            round(length(unique(tempo_hib_4dose$nodossier))/tempo_population*100,1),"\n",
            "intervalle de confiance : ",
            round(binom.test(dim(tempo_hib_4dose)[1],tempo_population)$conf.int*100,1)
      )
      return(sortie)
      
}
########################


#####################################################################
# fonction calcul CV 4 ans pour tempo_hib pour les cantons
f_hib4_canton <- function (df_hib){
      ############################################################
      # preparation
      # df_hib <- donnees[donnees$,] ; an <- 1997
      # restriction du df d'entree a l'annee voulue et au vaccin
      ddn_min <- as.Date("2007-01-01") # on veut les enfants ayants 24 mois 
      ddn_max <- as.Date("2011-12-31")
      tempo_hib <- df_hib[df_hib$datenaiss > ddn_min &
                          df_hib$datenaiss < ddn_max,] # temp_hib ne contient que les enfants de 24 mois l'annee choisie
      
      # calcul population denominateur avant restriction sur vaccins et age
      tempo_population <- length(unique(tempo_hib$nodossier))
      
      # on garde les vacci effectuees avant 4 ans (4 ans = 1460 jours) et en lien hib
      tempo_hib <- tempo_hib[tempo_hib$dateopv < (tempo_hib$datenaiss + 1460) &
                                   tempo_hib$vaccin_code %in% hib,]

      # on classe suivant le n° de dossier et l'anciennete de l'opv
      tempo_hib <- tempo_hib[order(tempo_hib$nodossier,tempo_hib$dateopv),]
      
      #############################################################
      # debut des calculs des nombres de vaccins
      
      # creation d'un df avec une seule repetition pour chaque dossier => nb de vacci 1 dose
      tempo_hib_1dose <- tempo_hib[!duplicated(tempo_hib$nodossier),]
      # pour 2 doses, on prend les autres vacci donc, sans "!", puis on réapplique
      # utilisation d'un df "pivot" tempo_hib_dose
      tempo_hib_dose <- tempo_hib[duplicated(tempo_hib$nodossier),] # toutes les vacci sauf la premiere dose
      tempo_hib_2dose <- tempo_hib_dose[!duplicated(tempo_hib_dose$nodossier),]
      # idem 3ieme dose
      tempo_hib_dose <- tempo_hib_dose[duplicated(tempo_hib_dose$nodossier),]
      tempo_hib_3dose <- tempo_hib_dose[!duplicated(tempo_hib_dose$nodossier),]
      # idem 4ieme dose
      tempo_hib_dose <- tempo_hib_dose[duplicated(tempo_hib_dose$nodossier),]
      tempo_hib_4dose <- tempo_hib_dose[!duplicated(tempo_hib_dose$nodossier),]
      
      sortie <- cat(
            "primovaccination hib dans le canton à 4 ans entre 2011 et 2015 en % : ",
            round(length(unique(tempo_hib_3dose$nodossier))/tempo_population*100,1),"\n",
            "intervalle de confiance : ",
            round(binom.test(dim(tempo_hib_3dose)[1],tempo_population)$conf.int*100,1),"\n",
            "rappel hib dans le canton à 4 ans entre 2011 et 2015 en % : ",
            round(length(unique(tempo_hib_4dose$nodossier))/tempo_population*100,1),"\n",
            "intervalle de confiance : ",
            round(binom.test(dim(tempo_hib_4dose)[1],tempo_population)$conf.int*100,1)
      )
      return(sortie)
      
}

######################################################################


for (c in levels(donnees$canton)){
      cat("
          ************************",c,"*********************************")
      
      tempo_canton <- donnees[donnees$canton == c,]
      cat("\n")
      f_hib24_canton(tempo_canton)
      cat("\n")
      f_hib4_canton(tempo_canton)
}


# ###### analyse des donnees
sink()

# moyenne <- read.csv2("sorties/hib/hib24_3_canton_pour moyenne.csv")
# 
# round(addmargins(moyenne,2,mean)[,6],1)
# 
# result <- matrix(round(as.numeric(addmargins(moyenne,2,mean)[,6]),1)[c(T,F,F)],2,6)
# 
# rownames(result) <- c("primovacci","rappel")
# colnames(result) <- levels(donnees$canton)
# write.csv2(result,"sorties/hib/result_cantons.csv")


sink("sorties/hib/hib_ana_primovacci.txt") # sorties vers le fichier specifie
cat("
      **********************************
      **********************************
      ****** ANALYSE AGE PRIMOVACCI ****
      **********************************
      **********************************")


      primovacci <- donnees[donnees$vaccin_code %in% hib,]

      # on classe suivant le n° de dossier et l'anciennete de l'opv
      primovacci <- primovacci[order(primovacci$nodossier,primovacci$dateopv),]
      
      #############################################################
      
      # contient "au moins 2 doses"
      primovacci <- primovacci[duplicated(primovacci$nodossier),] # toutes les vacci sauf la premiere dose
      # on reapplique, donc, "au moins 3 doses"
      primovacci <- primovacci[duplicated(primovacci$nodossier),]
      rappelvacci <- primovacci     # en attente, pour analyse age au rappel
      
      # mais pas plus, pour avoir l'age au dernier acte
      primovacci <- primovacci[!duplicated(primovacci$nodossier),]
      
      # analyse en tant que telle
      summary(primovacci$age_vacci_mois)
      
      # on reapplique, donc, "au moins  doses" c'est a dire rappel
      rappelvacci <- rappelvacci[duplicated(rappelvacci$nodossier),]
      rappelvacci <- rappelvacci[!duplicated(rappelvacci$nodossier),]

      # analyse de ce df
      png("sorties/hib/age-primovacci%01d.png", width=25 ,height=15, units="cm",res = 400)
      barplot(table(primovacci$age_vacci_mois[primovacci$age_vacci_mois < 48]),
           main = "Répartition de l'âge en mois à la primovaccination Hib,
pour les enfants nés entre le 1er janvier 1993
et le 31 décembre 2011, Corse-du-Sud",
           xlab = "âge en mois - coupure à 48 mois", ylab = "effectif") 
      
      barplot(table(rappelvacci$age_vacci_mois[rappelvacci$age_vacci_mois < 48]),
           main = "Répartition de l'âge en mois au rappel Hib,
pour les enfants nés entre le 1er janvier 1993
et le 31 décembre 2011, Corse-du-Sud",
           xlab = "âge en mois - coupure à 48 mois", ylab = "effectif") 
      dev.off()

sink()      


# fin analyse hib
####################################################################################################


#################################
#
# THAT'S ALL FOLKS !!!!!!!!!!!
#
#################################
