##################################
#
#   analyse donnees vacci CD 2A -- men
#
# V1  août 2017
##################################

rm(list = ls())

setwd("I:\\DIRECTION ACTION TERRITORIALE DE SANTE\\VSSE\\CIRE\\Dossiers MI\\Vaccination\\etude_CG2A\\analyse")
load("donnees/df_pour_analyse.RData")

men <- c("MEN C")


####################################################################
# analyse donnees 24 mois
###################################################################

donnees_CS24 <- read.csv2("donnees/donnees_CS_24_men.csv",header=F)

colnames(donnees_CS24) <- 2011:2016
rownames(donnees_CS24) <- c("national","Corse")


# analyse evolution men
for (i in 1:2) print(summary(lm(as.numeric(donnees_CS24[i,]) ~ as.numeric(2011:2016))))

# graphe evolution
couleurs <- c("red","blue")

png("sorties/men/evolution_CS24.png",width=25 ,height=15, units="cm",res = 400)
matplot(t(donnees_CS24),type="b",pch=1,xlab = "années", ylim = c(0,100), ylab = "taux de couverture vaccinal en %",
         main = "Evolution de la couverture vaccinale contre le méningocoque C à l'âge de 24 mois, 
         2011-2016 - résultats en %", xaxt="n",
        col = couleurs)
axis(1,at=1:6,labels=2011:2016)
abline(h = 95, col = "darkorange",lty = 2)
legend(5,14,legend = rownames(donnees_CS24), pch=1,col = couleurs, cex=0.8,bty="n",lwd=1)
legend(5,5,legend = "objectif de vaccination", pch="",col = "darkorange", cex=0.8,bty="n",lwd=1)
dev.off()
 

# fin analyse 24 mois men     
#################################

#####################################################################
# creation de la fonction calcul CV 24 mois pour tempo_men
f_men24 <- function (df_men,an){
      ############################################################
      # preparation
      # df_men <- donnees ; an <- 1997
      # restriction du df d'entree a l'annee voulue et au vaccin
      ddn_min <- as.Date(paste(an-2,"-01-01",sep="")) # on veut les enfants ayants 24 mois 
      ddn_max <- as.Date(paste(an-2,"-12-31",sep=""))
      tempo_men <- df_men[df_men$datenaiss > ddn_min &
                          df_men$datenaiss < ddn_max,] # temp_men ne contient que les enfants de 24 mois l'annee choisie
      
      tempo_population <- length(unique(tempo_men$nodossier)) # calcul population denominateur avant restriction sur vaccins et age a la vacci
      
      # on garde les vacci effectuees avant 24 mois (2 ans = 730 jours) et vacci men
      tempo_men <- tempo_men[tempo_men$dateopv < (tempo_men$datenaiss + 730)  &
                                   tempo_men$vaccin_code %in% men,]
      

      # on classe suivant le n° de dossier et l'anciennete de l'opv
      tempo_men <- tempo_men[order(tempo_men$nodossier,tempo_men$dateopv),]
      
      #############################################################
      # debut des calculs des nombres de vaccins
      
      # creation d'un df avec une seule repetition pour chaque dossier => nb de vacci 1 dose
      tempo_men_dose <- tempo_men[!duplicated(tempo_men$nodossier),]
      
      sortie <- list(
            round(dim(tempo_men_dose)[1] / tempo_population * 100,2),
            binom.test(dim(tempo_men_dose)[1],tempo_population)$conf.int,
            summary(as.numeric(tempo_men_dose$age_vacci_mois))
            )
      
      return(sortie)
      
}
########################

#####################################################################
# creation de la fonction calcul CV 4 ans pour tempo_men
f_men4 <- function (df_men,an){
      ############################################################
      # preparation
      
      # restriction du df d'entree a l'annee voulue et au vaccin
      ddn_min <- as.Date(paste(an-4,"-01-01",sep="")) # on veut les enfants ayants 4 ans l'annee choisie
      ddn_max <- as.Date(paste(an-4,"-12-31",sep=""))
      tempo_men <- df_men[df_men$datenaiss > ddn_min &
                          df_men$datenaiss < ddn_max,]
      
      tempo_population <- length(unique(tempo_men$nodossier)) # calcul population denominateur avant restriction sur vaccins
      
      # on garde les vacci effectuees avant 4 ans (= 1460 jours) et sur men
      tempo_men <- tempo_men[tempo_men$dateopv < (tempo_men$datenaiss + 1460) &
                                   tempo_men$vaccin_code %in% men,]
      
      
      # on classe suivant le n° de dossier et l'anciennete de l'opv
      tempo_men <- tempo_men[order(tempo_men$nodossier,tempo_men$dateopv),]
      
      #############################################################
      # debut des calculs des nombres de vaccines
      
      # creation d'un df avec une seule repetition pour chaque dossier => nb de vacci 1 dose
      tempo_men_dose <- tempo_men[!duplicated(tempo_men$nodossier),]
     
      sortie <- list(
            round(dim(tempo_men_dose)[1] / tempo_population * 100,2),
            binom.test(dim(tempo_men_dose)[1],tempo_population)$conf.int,
            summary(as.numeric(tempo_men_dose$age_vacci_mois))
            )
      
      return(sortie)
}
########################
###################################################
# boucle de sortie des resultats

sink("sorties/men/men_24.txt") # sorties vers le fichier specifie

print("###################################### 24 mois #########################################")

#remplissage des tableaux de sortie
sortie <- matrix(NA,7,4)
colnames(sortie) <- 2010:2013
rownames(sortie) <- c("CV en %","IC inf (%)","IC sup (%)","1er quartile","médiane","moyenne","3ième quartile")

for (i in 2010:2013) {
      tempo <- f_men24(donnees,i)
      sortie[1,i-2009] <- tempo[[1]] # CV
      sortie[2,i-2009] <- tempo[[2]][1]*100 ; sortie[3,i-2009] <- tempo[[2]][2]*100 # IC
      for (j in 4:7)
            sortie[j,i-2009] <- tempo[[3]][j-2] # repartition des ages en mois
}
print(round(sortie,2))
rm(tempo,sortie) # on efface les fichiers temporaires pour faire de la place

sink()

sink("sorties/men/men_4.txt") # sorties vers le fichier specifie
print("###################################### 4 ans #########################################")

#remplissage des tableaux de sortie
sortie <- matrix(NA,7,6)
colnames(sortie) <- 2010:2015
rownames(sortie) <- c("CV en %","IC inf (%)","IC sup (%)","1er quartile","médiane","moyenne","3ième quartile")

for (i in 2010:2015) {
      tempo <- f_men4(donnees,i)
      sortie[1,i-2009] <- tempo[[1]] # CV
      sortie[2,i-2009] <- tempo[[2]][1]*100 ; sortie[3,i-2009] <- tempo[[2]][2]*100 # IC
      for (j in 4:7)
            sortie[j,i-2009] <- tempo[[3]][j-2] # repartition des ages en mois
}
print(round(sortie,2))
rm(tempo,sortie) # on efface les fichiers temporaires pour faire de la place


sink()

sink("sorties/men/men_canton.txt") # sorties vers le fichier specifie
cat("
      ******************************
      ****** men *******************
      ******************************
      ****** PAR CANTONS ***********
      ******************************")

#####################################################################
# fonction calcul CV 24 mois pour tempo_men pour les cantons
f_men24_canton <- function (df_men){
      ############################################################
      # preparation
      # df_men <- donnees[donnees$,] ; an <- 1997
      # restriction du df d'entree a l'annee voulue et au vaccin
      ddn_min <- as.Date("2007-01-01") # on veut les enfants ayants 24 mois 
      ddn_max <- as.Date("2011-12-31")
      tempo_men <- df_men[df_men$datenaiss > ddn_min &
                          df_men$datenaiss < ddn_max,] # temp_men ne contient que les enfants de 24 mois l'annee choisie
      
      # calcul population denominateur avant restriction sur vaccins et age
      tempo_population <- length(unique(tempo_men$nodossier))
      
      # on garde les vacci effectuees avant 24 mois (2 ans = 730 jours) et en lien men
      tempo_men <- tempo_men[tempo_men$dateopv < (tempo_men$datenaiss + 730) &
                                   tempo_men$vaccin_code %in% men,]

      # on classe suivant le n° de dossier et l'anciennete de l'opv
      tempo_men <- tempo_men[order(tempo_men$nodossier,tempo_men$dateopv),]
      
      #############################################################
      # debut des calculs des nombres de vaccins
      
      # creation d'un df avec une seule repetition pour chaque dossier => nb de vacci 1 dose
      tempo_men_dose <- tempo_men[!duplicated(tempo_men$nodossier),]
      
      sortie <- cat(
      "primovaccination men dans le canton à 24 mois entre 2009 et 2013 en % : ",
      round(length(unique(tempo_men_dose$nodossier))/tempo_population*100,1),"\n",
      "intervalle de confiance : ",
      round(binom.test(dim(tempo_men_dose)[1],tempo_population)$conf.int*100,1))
      
      return(sortie)
      
}
########################


#####################################################################
# fonction calcul CV 4 ans pour tempo_men pour les cantons
f_men4_canton <- function (df_men){
      ############################################################
      # preparation
      # df_men <- donnees[donnees$,] ; an <- 1997
      # restriction du df d'entree a l'annee voulue et au vaccin
      ddn_min <- as.Date("2007-01-01") # on veut les enfants ayants 24 mois 
      ddn_max <- as.Date("2011-12-31")
      tempo_men <- df_men[df_men$datenaiss > ddn_min &
                          df_men$datenaiss < ddn_max,] # temp_men ne contient que les enfants de 24 mois l'annee choisie
      
      # calcul population denominateur avant restriction sur vaccins et age
      tempo_population <- length(unique(tempo_men$nodossier))
      
      # on garde les vacci effectuees avant 4 ans (4 ans = 1460 jours) et en lien men
      tempo_men <- tempo_men[tempo_men$dateopv < (tempo_men$datenaiss + 1460) &
                                   tempo_men$vaccin_code %in% men,]

      # on classe suivant le n° de dossier et l'anciennete de l'opv
      tempo_men <- tempo_men[order(tempo_men$nodossier,tempo_men$dateopv),]
      
      #############################################################
      # debut des calculs des nombres de vaccins
      
      # creation d'un df avec une seule repetition pour chaque dossier => nb de vacci 1 dose
      tempo_men_dose <- tempo_men[!duplicated(tempo_men$nodossier),]
      
      sortie <- cat(
            "primovaccination men dans le canton à 4 ans entre 2011 et 2015 en % : ",
            round(length(unique(tempo_men_dose$nodossier))/tempo_population*100,1),"\n",
            "intervalle de confiance : ",
            round(binom.test(dim(tempo_men_dose)[1],tempo_population)$conf.int*100,1)
      )
      return(sortie)
      
}

######################################################################


for (c in levels(donnees$canton)){
      cat("
          ************************",c,"*********************************")
      
      tempo_canton <- donnees[donnees$canton == c,]
      cat("\n")
      f_men24_canton(tempo_canton)
      cat("\n")
      f_men4_canton(tempo_canton)
}


# ###### analyse des donnees
sink()

# moyenne <- read.csv2("sorties/men/men24_3_canton_pour moyenne.csv")
# 
# round(addmargins(moyenne,2,mean)[,6],1)
# 
# result <- matrix(round(as.numeric(addmargins(moyenne,2,mean)[,6]),1)[c(T,F,F)],2,6)
# 
# rownames(result) <- c("primovacci","rappel")
# colnames(result) <- levels(donnees$canton)
# write.csv2(result,"sorties/men/result_cantons.csv")


sink("sorties/men/men_ana_primovacci.txt") # sorties vers le fichier specifie
cat("
      **********************************
      **********************************
      ****** ANALYSE AGE PRIMOVACCI ****
      **********************************
      **********************************")


      primovacci <- donnees[donnees$vaccin_code %in% men,]

      # on classe suivant le n° de dossier et l'anciennete de l'opv
      primovacci <- primovacci[order(primovacci$nodossier,primovacci$dateopv),]
      
      #############################################################
      
      # juste une dose
      primovacci <- primovacci[!duplicated(primovacci$nodossier),] 
      
      # analyse en tant que telle
      summary(primovacci$age_vacci_mois)
      # analyse de ce df
      png("sorties/men/age-primovacci%01d.png", width=30 ,height=15, units="cm",res = 400)
      barplot(table(primovacci$age_vacci_mois[primovacci$age_vacci_mois < 49]),
           main = "Répartition de l'âge en mois à la vaccination contre le méningocoque C,
pour les enfants nés entre le 1er janvier 1993
et le 31 décembre 2011, Corse-du-Sud",
           xlab = "âge en mois - coupure à 48 mois", ylab = "effectif",cex.axis = 0.8) 
      dev.off()

sink()      


# fin analyse men
####################################################################################################


#################################
#
# THAT'S ALL FOLKS !!!!!!!!!!!
#
#################################
