##################################
#
#   analyse donnees vacci CD 2A -- DTP
#
# V1  août 2017
##################################



# # nb de chaque vaccin
# table(donnees$vaccin_code[donnees$vaccin_code %in% dtp])
# 
# # annee vacci 
# plot(table(format(DTPolio_2015$dateopv,"%Y")))
# plot(table(DTPolio_2015$age_vacci)) 

load("donnees/df_pour_analyse.RData")

dtp <- c("D(A)TC(A)P","D(A)TP","DT","DTC","DTCP","DTCPHIB","DTCPHIBHB","DTP","DTP CA+HIB",
             "DTP+C(A)","DTTAB","P","T","TP")

sink("sorties/dtp/dtp.txt") # sorties vers le fichier specifie

print("***************************************************************************************")
print("*********************************** DTP ***********************************************")
print("***************************************************************************************")

# interet a separer les 3 valences 
tempo <- donnees[donnees$vaccin_code %in% dtp,] 

table(factor(tempo$vaccin_code)) # affichage de la repartition
print("pourcentage de vacci incomplètes DTP par rapport au total :")
paste(round(sum(table(factor(tempo$vaccin_code[tempo$vaccin_code %in% c("DT","DTC","P","T","TP")]))) /
      sum(table(factor(tempo$vaccin_code)))*100,3),"%")
png("sorties/dtp/repartition des vacci incompletes DTP.png")
plot(table(format(tempo$dateopv[tempo$vaccin_code %in% c("DT","DTC","P","T","TP")],"%Y")))
dev.off() # repartition graphique par annee des vacci incompletes

rm(tempo) # on fait de la place

# => aucun interet à separer, car pourcentage très faible, on restreint aux vacci avec les 3 valences
dtp <- c("D(A)TC(A)P","D(A)TP","DTCP","DTCPHIB","DTCPHIBHB","DTP","DTP CA+HIB", "DTP+C(A)")

####################################################################
# analyse donnees nationales
###################################################################

donnees_nat <- read.csv2("donnees/donnees_nationales_CS_24_DTP.csv",header=F)

colnames(donnees_nat) <- 1997:2015
rownames(donnees_nat) <- c("primovaccination DT", "rappel DT", "primovaccination P", "rappel P")


# analyse evolution DT
for (i in 1:4){
      print(rownames(donnees_nat)[i])
      print(summary(lm(as.numeric(donnees_nat[i,]) ~ as.numeric(1997:2015))))
      
}




# graphe evolution
couleurs <- c("blue","red","darkgreen", "black")

png("sorties/dtp/evolution_nat_CS24.png",width=25 ,height=15, units="cm",res = 200)
matplot(t(donnees_nat),type="b",pch=1:4,xlab = "années", ylim = c(85,100),ylab = "taux de couverture vaccinal en %",
         main = "Couverture vaccinale diphtérie, tétanos, poliomyélite à l'âge de 24 mois, 
         analyse des CS 24 - France, 1997-2015 - résultats en % - axe des ordonnées débutant à 85 %", xaxt="n",
        col = couleurs)
axis(1,at=1:19,labels=1997:2015)
legend(14,89,legend = rownames(donnees_nat), pch=1:4,col = couleurs, cex=0.8,bty="n",lwd=1)
abline(h = 95, col = "darkorange",lty = 2)
dev.off()
 

# fin analyse nationale DTP     
#################################

#####################################################################
# creation de la fonction calcul CV 24 mois pour tempo_dtp
f_dtp24 <- function (df_dtp,an){
      ############################################################
      # preparation
      # df_dtp <- donnees ; an <- 1997
      # restriction du df d'entree a l'annee voulue et au vaccin
      ddn_min <- as.Date(paste(an-2,"-01-01",sep="")) # on veut les enfants ayants 24 mois 
      ddn_max <- as.Date(paste(an-2,"-12-31",sep=""))
      tempo_dtp <- df_dtp[df_dtp$datenaiss > ddn_min &
                          df_dtp$datenaiss < ddn_max,] # temp_dtp ne contient que les enfants de 24 mois l'annee choisie
      
      tempo_population <- length(unique(tempo_dtp$nodossier)) # calcul population denominateur avant restriction sur vaccins et age a la vacci
      
      # on garde les vacci effectuees avant 24 mois (2 ans = 730 jours) et vacci DTP
      tempo_dtp <- tempo_dtp[tempo_dtp$dateopv < (tempo_dtp$datenaiss + 730)  &
                                   tempo_dtp$vaccin_code %in% dtp,]
      

      # on classe suivant le n° de dossier et l'anciennete de l'opv
      tempo_dtp <- tempo_dtp[order(tempo_dtp$nodossier,tempo_dtp$dateopv),]
      
      #############################################################
      # debut des calculs des nombres de vaccins
      
      # creation d'un df avec une seule repetition pour chaque dossier => nb de vacci 1 dose
      tempo_dtp_1dose <- tempo_dtp[!duplicated(tempo_dtp$nodossier),]
      # pour 2 doses, on prend les autres vacci donc, sans "!", puis on réapplique
      # utilisation d'un df "pivot" tempo_dtp_dose
      tempo_dtp_dose <- tempo_dtp[duplicated(tempo_dtp$nodossier),] # toutes les vacci sauf la premiere dose
      tempo_dtp_2dose <- tempo_dtp_dose[!duplicated(tempo_dtp_dose$nodossier),]
      # idem 3ieme dose
      tempo_dtp_dose <- tempo_dtp_dose[duplicated(tempo_dtp_dose$nodossier),]
      tempo_dtp_3dose <- tempo_dtp_dose[!duplicated(tempo_dtp_dose$nodossier),]
      # idem 4ieme dose
      tempo_dtp_dose <- tempo_dtp_dose[duplicated(tempo_dtp_dose$nodossier),]
      tempo_dtp_4dose <- tempo_dtp_dose[!duplicated(tempo_dtp_dose$nodossier),]
      
      
      
      sortie1 <- list(
            round(dim(tempo_dtp_1dose)[1] / tempo_population * 100,2),
            binom.test(dim(tempo_dtp_1dose)[1],tempo_population)$conf.int,
            summary(as.numeric(tempo_dtp_1dose$age_vacci_mois[tempo_dtp_1dose$acte_code == "P1"]))
            )
      sortie2 <- list(
            round(dim(tempo_dtp_2dose)[1] / tempo_population * 100,2),
            binom.test(dim(tempo_dtp_2dose)[1],tempo_population)$conf.int,
            summary(as.numeric(tempo_dtp_2dose$age_vacci_mois[tempo_dtp_2dose$acte_code == "P2"]))
            )
      sortie3 <- list(
            round(dim(tempo_dtp_3dose)[1] / tempo_population * 100,2),
            binom.test(dim(tempo_dtp_3dose)[1],tempo_population)$conf.int,
            summary(as.numeric(tempo_dtp_3dose$age_vacci_mois[tempo_dtp_3dose$acte_code == "P3"]))
            )
      sortie4 <- list(
            round(dim(tempo_dtp_4dose)[1] / tempo_population * 100,2),
            binom.test(dim(tempo_dtp_4dose)[1],tempo_population)$conf.int,
            summary(as.numeric(tempo_dtp_4dose$age_vacci_mois[tempo_dtp_4dose$acte_code == "R01"]))
            )
      sortie <- list(sortie1,sortie2,sortie3,sortie4) # liste de listes ...
      return(sortie)
      
}
########################

#####################################################################
# creation de la fonction calcul CV 4 ans pour tempo_dtp
f_dtp4 <- function (df_dtp,an){
      ############################################################
      # preparation
      
      # restriction du df d'entree a l'annee voulue et au vaccin
      ddn_min <- as.Date(paste(an-4,"-01-01",sep="")) # on veut les enfants ayants 4 ans l'annee choisie
      ddn_max <- as.Date(paste(an-4,"-12-31",sep=""))
      tempo_dtp <- df_dtp[df_dtp$datenaiss > ddn_min &
                          df_dtp$datenaiss < ddn_max,]
      
      tempo_population <- length(unique(tempo_dtp$nodossier)) # calcul population denominateur avant restriction sur vaccins
      
      # on garde les vacci effectuees avant 4 ans (= 1460 jours) et sur dtp
      tempo_dtp <- tempo_dtp[tempo_dtp$dateopv < (tempo_dtp$datenaiss + 1460) &
                                   tempo_dtp$vaccin_code %in% dtp,]
      
      
      # on classe suivant le n° de dossier et l'anciennete de l'opv
      tempo_dtp <- tempo_dtp[order(tempo_dtp$nodossier,tempo_dtp$dateopv),]
      
      #############################################################
      # debut des calculs des nombres de vaccines
      
      # creation d'un df avec une seule repetition pour chaque dossier => nb de vacci 1 dose
      tempo_dtp_1dose <- tempo_dtp[!duplicated(tempo_dtp$nodossier),]
      # pour 2 doses, on prend les autres vacci donc, sans "!", puis on réapplique
      # utilisation d'un df "pivot" tempo_dtp_dose
      tempo_dtp_dose <- tempo_dtp[duplicated(tempo_dtp$nodossier),] # toutes les vacci sauf la premiere dose
      tempo_dtp_2dose <- tempo_dtp_dose[!duplicated(tempo_dtp_dose$nodossier),]
      # idem 3ieme dose
      tempo_dtp_dose <- tempo_dtp_dose[duplicated(tempo_dtp_dose$nodossier),]
      tempo_dtp_3dose <- tempo_dtp_dose[!duplicated(tempo_dtp_dose$nodossier),]
      # idem 4ieme dose
      tempo_dtp_dose <- tempo_dtp_dose[duplicated(tempo_dtp_dose$nodossier),]
      tempo_dtp_4dose <- tempo_dtp_dose[!duplicated(tempo_dtp_dose$nodossier),]
      
      sortie1 <- list(
            round(dim(tempo_dtp_1dose)[1] / tempo_population * 100,2),
            binom.test(dim(tempo_dtp_1dose)[1],tempo_population)$conf.int,
            summary(as.numeric(tempo_dtp_1dose$age_vacci_mois[tempo_dtp_1dose$acte_code == "P1"]))
            )
      sortie2 <- list(
            round(dim(tempo_dtp_2dose)[1] / tempo_population * 100,2),
            binom.test(dim(tempo_dtp_2dose)[1],tempo_population)$conf.int,
            summary(as.numeric(tempo_dtp_2dose$age_vacci_mois[tempo_dtp_2dose$acte_code == "P2"]))
            )
      sortie3 <- list(
            round(dim(tempo_dtp_3dose)[1] / tempo_population * 100,2),
            binom.test(dim(tempo_dtp_3dose)[1],tempo_population)$conf.int,
            summary(as.numeric(tempo_dtp_3dose$age_vacci_mois[tempo_dtp_3dose$acte_code == "P3"]))
            )
      sortie4 <- list(
            round(dim(tempo_dtp_4dose)[1] / tempo_population * 100,2),
            binom.test(dim(tempo_dtp_4dose)[1],tempo_population)$conf.int,
            summary(as.numeric(tempo_dtp_4dose$age_vacci_mois[tempo_dtp_4dose$acte_code == "R01"]))
            )
      sortie <- list(sortie1,sortie2,sortie3,sortie4) # liste de listes ...
      return(sortie)
}
########################
###################################################
# boucle de sortie des resultats



print("###################################### 24 mois #########################################")

#remplissage des tableaux de sortie
for (d in 3:4){ # d pour dose
      print(paste("---------------------------",d," dose(s) ----------------------------"))
      
      sortie <- matrix(NA,7,19)
      colnames(sortie) <- 1995:2013
      rownames(sortie) <- c("CV en %","IC inf (%)","IC sup (%)","1er quartile","médiane","moyenne","3ième quartile")
   
      for (i in 1995:2013) {
            tempo <- f_dtp24(donnees,i)[[d]] # on ne recupere que la "d"ieme liste dans la sortie
            sortie[1,i-1994] <- tempo[[1]] # CV
            sortie[2,i-1994] <- tempo[[2]][1]*100 ; sortie[3,i-1994] <- tempo[[2]][2]*100 # IC
            for (j in 4:7)
                  sortie[j,i-1994] <- tempo[[3]][j-2] # repartition des ages en mois
      }
      print(round(sortie,2))
      rm(tempo,sortie) # on efface les fichiers temporaires pour faire de la place
}




print("###################################### 4 ans #########################################")

#remplissage des tableaux de sortie
for (d in 3:4){ # d pour dose
      print(paste("---------------------------",d," dose(s) ----------------------------"))
      
      sortie <- matrix(NA,7,19)
      colnames(sortie) <- 1997:2015
      rownames(sortie) <- c("CV en %","IC inf (%)","IC sup (%)","1er quartile","médiane","moyenne","3ième quartile")
   
      for (i in 1997:2015) {
            tempo <- f_dtp4(donnees,i)[[d]] # on ne recupere que la "d"ieme liste dans la sortie
            sortie[1,i-1996] <- tempo[[1]] # CV
            sortie[2,i-1996] <- tempo[[2]][1]*100 ; sortie[3,i-1996] <- tempo[[2]][2]*100 # IC
            for (j in 4:7)
                  sortie[j,i-1996] <- tempo[[3]][j-2] # repartition des ages en mois
      }
      print(round(sortie,2))
      rm(tempo,sortie) # on efface les fichiers temporaires pour faire de la place
}




cat("
      ******************************
      ****** DTP *******************
      ******************************
      ****** PAR CANTONS ***********
      ******************************")

#####################################################################
# fonction calcul CV 24 mois pour tempo_dtp pour les cantons
f_dtp24_canton <- function (df_dtp){
      ############################################################
      # preparation
      # df_dtp <- donnees[donnees$,] ; an <- 1997
      # restriction du df d'entree a l'annee voulue et au vaccin
      ddn_min <- as.Date("2007-01-01") # on veut les enfants ayants 24 mois 
      ddn_max <- as.Date("2011-12-31")
      tempo_dtp <- df_dtp[df_dtp$datenaiss > ddn_min &
                          df_dtp$datenaiss < ddn_max,] # temp_dtp ne contient que les enfants de 24 mois l'annee choisie
      
      # calcul population denominateur avant restriction sur vaccins et age
      tempo_population <- length(unique(tempo_dtp$nodossier))
      
      # on garde les vacci effectuees avant 24 mois (2 ans = 730 jours) et en lien dtp
      tempo_dtp <- tempo_dtp[tempo_dtp$dateopv < (tempo_dtp$datenaiss + 730) &
                                   tempo_dtp$vaccin_code %in% dtp,]

      # on classe suivant le n° de dossier et l'anciennete de l'opv
      tempo_dtp <- tempo_dtp[order(tempo_dtp$nodossier,tempo_dtp$dateopv),]
      
      #############################################################
      # debut des calculs des nombres de vaccins
      
      # creation d'un df avec une seule repetition pour chaque dossier => nb de vacci 1 dose
      tempo_dtp_1dose <- tempo_dtp[!duplicated(tempo_dtp$nodossier),]
      # pour 2 doses, on prend les autres vacci donc, sans "!", puis on réapplique
      # utilisation d'un df "pivot" tempo_dtp_dose
      tempo_dtp_dose <- tempo_dtp[duplicated(tempo_dtp$nodossier),] # toutes les vacci sauf la premiere dose
      tempo_dtp_2dose <- tempo_dtp_dose[!duplicated(tempo_dtp_dose$nodossier),]
      # idem 3ieme dose
      tempo_dtp_dose <- tempo_dtp_dose[duplicated(tempo_dtp_dose$nodossier),]
      tempo_dtp_3dose <- tempo_dtp_dose[!duplicated(tempo_dtp_dose$nodossier),]
      # idem 4ieme dose
      tempo_dtp_dose <- tempo_dtp_dose[duplicated(tempo_dtp_dose$nodossier),]
      tempo_dtp_4dose <- tempo_dtp_dose[!duplicated(tempo_dtp_dose$nodossier),]
      
      sortie <- cat(
            "primovaccination DTP dans le canton à 24 mois entre 2009 et 2013 en % : ",
            round(length(unique(tempo_dtp_3dose$nodossier))/tempo_population*100,1),"\n",
            "intervalle de confiance : ",
            round(binom.test(dim(tempo_dtp_3dose)[1],tempo_population)$conf.int*100,1),"\n",
            "rappel DTP dans le canton à 24 mois entre 2009 et 2013 en % : ",
            round(length(unique(tempo_dtp_4dose$nodossier))/tempo_population*100,1),"\n",
            "intervalle de confiance : ",
            round(binom.test(dim(tempo_dtp_4dose)[1],tempo_population)$conf.int*100,1)
      )
      return(sortie)
      
}
########################


#####################################################################
# fonction calcul CV 4 ans pour tempo_dtp pour les cantons
f_dtp4_canton <- function (df_dtp){
      ############################################################
      # preparation
      # df_dtp <- donnees[donnees$,] ; an <- 1997
      # restriction du df d'entree a l'annee voulue et au vaccin
      ddn_min <- as.Date("2007-01-01") # on veut les enfants ayants 24 mois 
      ddn_max <- as.Date("2011-12-31")
      tempo_dtp <- df_dtp[df_dtp$datenaiss > ddn_min &
                          df_dtp$datenaiss < ddn_max,] # temp_dtp ne contient que les enfants de 24 mois l'annee choisie
      
      # calcul population denominateur avant restriction sur vaccins et age
      tempo_population <- length(unique(tempo_dtp$nodossier))
      
      # on garde les vacci effectuees avant 4 ans (4 ans = 1460 jours) et en lien dtp
      tempo_dtp <- tempo_dtp[tempo_dtp$dateopv < (tempo_dtp$datenaiss + 1460) &
                                   tempo_dtp$vaccin_code %in% dtp,]

      # on classe suivant le n° de dossier et l'anciennete de l'opv
      tempo_dtp <- tempo_dtp[order(tempo_dtp$nodossier,tempo_dtp$dateopv),]
      
      #############################################################
      # debut des calculs des nombres de vaccins
      
      # creation d'un df avec une seule repetition pour chaque dossier => nb de vacci 1 dose
      tempo_dtp_1dose <- tempo_dtp[!duplicated(tempo_dtp$nodossier),]
      # pour 2 doses, on prend les autres vacci donc, sans "!", puis on réapplique
      # utilisation d'un df "pivot" tempo_dtp_dose
      tempo_dtp_dose <- tempo_dtp[duplicated(tempo_dtp$nodossier),] # toutes les vacci sauf la premiere dose
      tempo_dtp_2dose <- tempo_dtp_dose[!duplicated(tempo_dtp_dose$nodossier),]
      # idem 3ieme dose
      tempo_dtp_dose <- tempo_dtp_dose[duplicated(tempo_dtp_dose$nodossier),]
      tempo_dtp_3dose <- tempo_dtp_dose[!duplicated(tempo_dtp_dose$nodossier),]
      # idem 4ieme dose
      tempo_dtp_dose <- tempo_dtp_dose[duplicated(tempo_dtp_dose$nodossier),]
      tempo_dtp_4dose <- tempo_dtp_dose[!duplicated(tempo_dtp_dose$nodossier),]
      
      sortie <- cat(
            "primovaccination DTP dans le canton à 4 ans entre 2011 et 2015 en % : ",
            round(length(unique(tempo_dtp_3dose$nodossier))/tempo_population*100,1),"\n",
            "intervalle de confiance : ",
            round(binom.test(dim(tempo_dtp_3dose)[1],tempo_population)$conf.int*100,1),"\n",
            "rappel DTP dans le canton à 4 ans entre 2011 et 2015 en % : ",
            round(length(unique(tempo_dtp_4dose$nodossier))/tempo_population*100,1),"\n",
            "intervalle de confiance : ",
            round(binom.test(dim(tempo_dtp_4dose)[1],tempo_population)$conf.int*100,1)
      )
      return(sortie)
      
}

######################################################################


for (c in levels(donnees$canton)){
      cat("
          ************************",c,"*********************************")
      
      tempo_canton <- donnees[donnees$canton == c,]
      cat("\n")
      f_dtp24_canton(tempo_canton)
      cat("\n")
      f_dtp4_canton(tempo_canton)
}


# ###### analyse des donnees


moyenne <- read.csv2("sorties/dtp/dtp24_3_canton_pour moyenne.csv")

round(addmargins(moyenne,2,mean)[,6],1)

result <- matrix(round(as.numeric(addmargins(moyenne,2,mean)[,6]),1)[c(T,F,F)],2,6)

rownames(result) <- c("primovacci","rappel")
colnames(result) <- levels(donnees$canton)
write.csv2(result,"sorties/dtp/result_cantons.csv")



cat("
      **********************************
      **********************************
      ****** ANALYSE AGE PRIMOVACCI ****
      **********************************
      **********************************")


      primovacci <- donnees[donnees$vaccin_code %in% dtp,]

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
      # analyse de ce df
      png("sorties/dtp/age-primovacci%01d.png")
      plot(table(primovacci$age_vacci_mois[primovacci$age_vacci_mois < 36]),
           main = "Répartition de l'âge en mois à la primovaccination DTP,
pour les enfants nés entre le 1er janvier 1993
et le 31 décembre 2011, Corse-du-Sud",
           xlab = "âge en mois - coupure à 36 mois", ylab = "effectif") 
      dev.off()

      

cat("
      **********************************
      **********************************
      ******** ANALYSE AGE RAPPEL ******
      **********************************
      **********************************")


      #############################################################
      
      # on reapplique une derniere fois pour enlever la derniere injection de la primovacci
      rappelvacci <- rappelvacci[duplicated(rappelvacci$nodossier),] 
      
      # mais pas plus, pour avoir l'age au dernier acte
      rappelvacci <- rappelvacci[!duplicated(rappelvacci$nodossier),]
      
      # analyse en tant que telle
      summary(rappelvacci$age_vacci_mois)
      # analyse de ce df
      png("sorties/dtp/age-rappelvacci%01d.png")
      plot(table(rappelvacci$age_vacci_mois[rappelvacci$age_vacci_mois < 60]),
           main = "Répartition de l'âge en mois lors du rappel DTP,
pour les enfants nés entre le 1er janvier 1993
et le 31 décembre 2011, Corse-du-Sud",
           xlab = "âge en mois - coupure à 60 mois", ylab = "effectif") 
      dev.off()

sink()









# fin analyse DTPolio
####################################################################################################


#################################
#
# THAT'S ALL FOLKS !!!!!!!!!!!
#
#################################
