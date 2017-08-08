##################################
#
#   analyse donnees vacci CD 2A -- script d'analyse
#
# V0.1  avril 2017
##################################

rm(list=ls())

setwd("I:\\DIRECTION ACTION TERRITORIALE DE SANTE\\VSSE\\CIRE\\Dossiers MI\\Vaccination\\etude_CG2A\\analyse")
# import des donnees
load("donnees/df_pour_analyse.RData")

library(ggplot2)
library(reshape2)



######################################################################
# donnees population 2A
######################################################################

pop2A <- read.csv2("donnees/pop2A.csv")

# aggregation des donnees
pop2A_2013 <- aggregate(
      pop2A$population,
      list(pop2A$age),
      sum
)
colnames(pop2A_2013) <- c("age","population")


# pop2A <- read.csv2("donnees/pourGuillaumeIVS.csv")
# pop2A_2010 <- pop2A[pop2A$annee_recens == 2010,1:2]
# pop2A_2011 <- pop2A[pop2A$annee_recens == 2011,1:2]
# pop2A_2012 <- pop2A[pop2A$annee_recens == 2012,1:2]

# fin population 2A
#########################################################################

#########################################################################
# creation d'un df cantons
#########################################################################

# aggregation des donnees
pop2A_2013_cantons <- aggregate(
      pop2A$population,
      list(pop2A$age,pop2A$com),
      sum
)
colnames(pop2A_2013_cantons) <- c("age","commune","population")

communes_cantons <- read.csv2("donnees/communes_cantons.csv")
communes_cantons$Code.Insee <- as.numeric(substr(communes_cantons$Code.Insee,3,5))

# on fusionne les 2 et on ne garde que les colonnes utiles
pop2A_2013_cantons <- merge(pop2A_2013_cantons,communes_cantons,by.x = "commune",by.y = "Code.Insee")
pop2A_2013_cantons <- pop2A_2013_cantons[,c(1:5,8)]

# on somme par canton et par âge
cantons_2013 <- aggregate(pop2A_2013_cantons$population,list(pop2A_2013_cantons$age,pop2A_2013_cantons$canton),sum)
colnames(cantons_2013) <- c("age","cantons","population")

jeunes <- cantons_2013[cantons_2013$age %in% c(3,4),]
aggregate(jeunes$population,list(jeunes$cantons),sum)



# fin creation d'un df cantons
#########################################################################



##########################################################
# regroupement de vaccins
##########################################################
menC <- c("MEN C")
menACYW <- c("MEN C","M ACYW","M ACYW135") 

coq <- c("C","D(A)TC(A)P","DTC","DTCP","DTCPHIB","DTCPHIBHB","DTP CA+HIB","DTP+C(A)")
hib <- c("DTCPHIB","DTCPHIBHB","DTP CA+HIB","HIB")
vhb <- c("DTCPHIBHB","HA+HB","HB")
ror <- c("O","ROR","ROU","RRUB","RUB")


# fin regroupement codes vaccins
###########################################################################

###########################################################################
# activité brute du CV
###########################################################################

summary(donnees$dateopv)
plot(table(donnees$dateopv))

# p <- ggplot(donnees)
# p + geom_bar(aes(x=dateopv,colour = "blue"))
# p + geom_histogram(aes(x=dateopv,colour = "blue"))
# 
# donnees$age_fin_2016 <- as.numeric(floor((as.Date("2016-12-31")-donnees$datenaiss)/365.25))
# png("sorties/age_fin_2016.png",width=25 ,height=15, units="cm",res = 200)
# p + geom_bar(aes(x=age_fin_2016))
# dev.off()

# fin activite brute du centre de vacci
############################################################################

############################################################################
# analyse de la population
############################################################################
# creation d'un df nomme pop : population dans la base sans double compte (ie une ligne un vaccin)
pop <- donnees[,c(1,2,3,4,5,30)]
pop <- pop[order(pop$nodossier),]
pop <- unique(pop) 
summary(pop$datenaiss)

png("sorties/ddn_unique.png",width=25 ,height=15, units="cm",res = 200)
ddn <- ggplot(pop)
ddn + geom_bar(aes(x=datenaiss))
dev.off()

plot(table(pop$datenaiss))

# evolution du taux de completude de la base suivant les ages de 0 a 50 ans

evolution <- matrix(0,1,11) # on initialise la matrice de sortie avec des 0 et une ligne
colnames(evolution) <- 0:10
rownames(evolution) <- "completude en %"
for (ag in 0:10){
      evolution[ag+1] <- comptage(donnees,2015,ag)
      evolution[ag+1] <- round(evolution[ag+1] / pop2A_2013$population[pop2A_2013$age == ag] * 100,1)
}


plot(t(evolution))
compl <- as.data.frame(t(evolution))
compl <- as.data.frame(compl[3:26,])
summary(compl)

("sexe-ratio H/F")
table(pop$sexe)[3]/table(pop$sexe)[1]

("répartition par canton")
table(pop$canton)/dim(pop[!is.na(pop$canton),])[1]*100

summary(pop$canton)


 # fin analyse population
####################################################################################################





####################################################################################################
# ANALYSE COQUELUCHE
#################################################################################
table(donnees$vaccin_code)
summary(donnees$dateopv[donnees$vaccin_code == "C"])
plot(table(format(donnees$dateopv[donnees$vaccin_code == "C"],"%Y")))
################
##########################################################################################


# creation de la fonction calcul CV pour tempo_coq
f_coq <- function (df_coq,an){
      ############################################################
      # preparation
      
      # restriction du df d'entree a l'annee voulue et au vaccin
      ddn_min <- as.Date(paste(an-3,"-01-01",sep="")) # on veut les enfants ayants 3/4 ans l'annee choisie
      ddn_max <- as.Date(paste(an-3,"-12-31",sep=""))
      tempo_coq <- df_coq[df_coq$datenaiss > ddn_min &
                          df_coq$datenaiss < ddn_max,]
      
      tempo_population <- comptage(tempo_coq,an,3) # calcul population denominateur avant restriction sur vaccines
      tempo_coq <- tempo_coq[tempo_coq$vaccin_code %in% coq,]
      
      # on classe suivant le n° de dossier et l'anciennete de l'opv
      tempo_coq <- tempo_coq[order(tempo_coq$nodossier,tempo_coq$dateopv),]
      
      #############################################################
      # debut des calculs des nombres de vaccins
      
      # creation d'un df avec une seule repetition pour chaque dossier => nb de vacci 1 dose
      tempo_coq_1dose <- tempo_coq[!duplicated(tempo_coq$nodossier),]
      # pour 2 doses, on prend les autres vacci donc, sans "!", puis on réapplique
      # utilisation d'un df "pivot" tempo_coq_dose
      tempo_coq_dose <- tempo_coq[duplicated(tempo_coq$nodossier),] # toutes les vacci sauf la premiere dose
      tempo_coq_2dose <- tempo_coq_dose[!duplicated(tempo_coq_dose$nodossier),]
      # idem 3ieme dose
      tempo_coq_dose <- tempo_coq_dose[duplicated(tempo_coq_dose$nodossier),]
      tempo_coq_3dose <- tempo_coq_dose[!duplicated(tempo_coq_dose$nodossier),]
      # idem 4ieme dose
      tempo_coq_dose <- tempo_coq_dose[duplicated(tempo_coq_dose$nodossier),]
      tempo_coq_4dose <- tempo_coq_dose[!duplicated(tempo_coq_dose$nodossier),]
      
      sortie1 <- list(
            round(dim(tempo_coq_1dose)[1] / tempo_population * 100,2),
            binom.test(dim(tempo_coq_1dose)[1],tempo_population)$conf.int,
            summary(as.numeric(tempo_coq_1dose$age_vacci_mois[tempo_coq_1dose$acte_code == "P1"]))
            )
      sortie2 <- list(
            round(dim(tempo_coq_2dose)[1] / tempo_population * 100,2),
            binom.test(dim(tempo_coq_2dose)[1],tempo_population)$conf.int,
            summary(as.numeric(tempo_coq_2dose$age_vacci_mois[tempo_coq_2dose$acte_code == "P2"]))
            )
      sortie3 <- list(
            round(dim(tempo_coq_3dose)[1] / tempo_population * 100,2),
            binom.test(dim(tempo_coq_3dose)[1],tempo_population)$conf.int,
            summary(as.numeric(tempo_coq_3dose$age_vacci_mois[tempo_coq_3dose$acte_code == "P3"]))
            )
      sortie4 <- list(
            round(dim(tempo_coq_4dose)[1] / tempo_population * 100,2),
            binom.test(dim(tempo_coq_4dose)[1],tempo_population)$conf.int,
            summary(as.numeric(tempo_coq_4dose$age_vacci_mois[tempo_coq_4dose$acte_code == "R01"]))
            )
      sortie <- list(sortie1,sortie2,sortie3,sortie4) # liste de listes ...
      return(sortie)
}


###################################################
# boucle de sortie des resultats

print("*****************************
      ****** COQ *****************
      ******************************")

#remplissage des tableaux de sortie
for (d in 1:4){ # d pour dose
      print(paste("---------------------------",d," dose(s) ----------------------------"))
      
      sortie <- matrix(NA,7,6)
      colnames(sortie) <- 2010:2015
      rownames(sortie) <- c("CV en %","IC inf (%)","IC sup (%)","1er quartile","médiane","moyenne","3ième quartile")
   
      for (i in 2010:2015) {
            tempo <- f_coq(donnees,i)[[d]] # on ne recupere que la "d"ieme liste dans la sortie
            sortie[1,i-2009] <- tempo[[1]] # CV
            sortie[2,i-2009] <- tempo[[2]][1]*100 ; sortie[3,i-2009] <- tempo[[2]][2]*100 # IC
            for (j in 4:7)
                  sortie[j,i-2009] <- tempo[[3]][j-2] # repartition des ages en mois
      }
      print(round(sortie,2))
      rm(tempo,sortie) # on efface les fichiers temporaires pour faire de la place
}





print("*****************************
      ****** COQ *****************
      ******************************
      ****** PAR CANTONS ***********
      ******************************")



for (c in levels(donnees$canton)){
      print(paste("*****************************",c,"************************************************************"))
      tempo_canton <- donnees[donnees$canton == c,]
      for (d in 1:4) {
            print(paste("---------------------------",d," dose(s) ----------------------------"))
            sortie <- matrix(NA,7,6)
            colnames(sortie) <- 2010:2015
            rownames(sortie) <- c("CV en %","IC inf (%)","IC sup (%)","1er quartile","médiane","moyenne","3ième quartile")
      
            #remplissage du tableau de sortie
            for (i in 2010:2015) {
                  tempo <- f_coq(tempo_canton,i)[[d]] # on ne recupere que la "d"ieme liste dans la sortie
                  sortie[1,i-2009] <- tempo[[1]] # CV
                  sortie[2,i-2009] <- tempo[[2]][1]*100 ; sortie[3,i-2009] <- tempo[[2]][2]*100 # IC
                  for (j in 4:7)
                        sortie[j,i-2009] <- tempo[[3]][j-2] # repartition des ages en mois
            }
            print(round(sortie,2))
            # rm(tempo,tempo_canton,sortie) # on efface les fichiers temporaires pour faire de la place
      }
}



# fin analyse coqueluche
####################################################################################################

####################################################################################################
# ANALYSE HAEMOPHILUS INFLUENZA DE TYPE B
####################################################################################################

# creation de la fonction calcul CV pour tempo_hib
f_hib <- function (df_hib,an){
      ############################################################
      # preparation
      
      # restriction du df d'entree a l'annee voulue et au vaccin
      ddn_min <- as.Date(paste(an-3,"-01-01",sep="")) # on veut les enfants ayants 3/4 ans l'annee choisie
      ddn_max <- as.Date(paste(an-3,"-12-31",sep=""))
      tempo_hib <- df_hib[df_hib$datenaiss > ddn_min &
                          df_hib$datenaiss < ddn_max,]
      
      tempo_population <- comptage(tempo_hib,an,3) # calcul population denominateur avant restriction sur vaccines
      tempo_hib <- tempo_hib[tempo_hib$vaccin_code %in% hib,]
      
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
      
      sortie <- list(
            ("hib, 1iere dose, resume des ages"), summary(as.numeric(tempo_hib_1dose$age_vacci_mois[tempo_hib_1dose$acte_code == "P1"])),
            ("hib, 2ieme dose, resume des ages"), summary(as.numeric(tempo_hib_2dose$age_vacci_mois[tempo_hib_2dose$acte_code == "P2"])),
            ("hib, 3ieme dose, resume des ages"), summary(as.numeric(tempo_hib_3dose$age_vacci_mois[tempo_hib_3dose$acte_code == "P3"])),
            ("hib, 4ieme dose, resume des ages"), summary(as.numeric(tempo_hib_4dose$age_vacci_mois[tempo_hib_4dose$acte_code == "R01"])),
            ("hib, CV, 1 dose"), dim(tempo_hib_1dose)[1] / tempo_population * 100,
            ("hib, CV, 2 doses"), dim(tempo_hib_2dose)[1] / tempo_population * 100,
            ("hib, CV, 3 doses"), dim(tempo_hib_3dose)[1] / tempo_population * 100,
            ("hib, CV, 4 doses"), dim(tempo_hib_4dose)[1] / tempo_population * 100
     )
      return(sortie)
}

for (i in 2010:2015) {
      print(paste("*************",i,"***************"))
      print(f_hib(donnees,i))
}


# fin analyse Haemophilus influenza de type b
####################################################################################################

####################################################################################################
# ANALYSE HEPATITE B
####################################################################################################

# creation de la fonction calcul CV pour tempo_vhb
f_vhb <- function (df_vhb,an){
      ############################################################
      # preparation
      
      # restriction du df d'entree a l'annee voulue et au vaccin
      ddn_min <- as.Date(paste(an-3,"-01-01",sep="")) # on veut les enfants ayants 3/4 ans l'annee choisie
      ddn_max <- as.Date(paste(an-3,"-12-31",sep=""))
      tempo_vhb <- df_vhb[df_vhb$datenaiss > ddn_min &
                          df_vhb$datenaiss < ddn_max,]
      
      tempo_population <- comptage(tempo_vhb,an,3) # calcul population denominateur avant restriction sur vaccins
      tempo_vhb <- tempo_vhb[tempo_vhb$vaccin_code %in% vhb,]
      
      # on classe suivant le n° de dossier et l'anciennete de l'opv
      tempo_vhb <- tempo_vhb[order(tempo_vhb$nodossier,tempo_vhb$dateopv),]
      
      #############################################################
      # debut des calculs des nombres de vaccines
      
      # creation d'un df avec une seule repetition pour chaque dossier => nb de vacci 1 dose
      tempo_vhb_1dose <- tempo_vhb[!duplicated(tempo_vhb$nodossier),]
      # pour 2 doses, on prend les autres vacci donc, sans "!", puis on réapplique
      # utilisation d'un df "pivot" tempo_vhb_dose
      tempo_vhb_dose <- tempo_vhb[duplicated(tempo_vhb$nodossier),] # toutes les vacci sauf la premiere dose
      tempo_vhb_2dose <- tempo_vhb_dose[!duplicated(tempo_vhb_dose$nodossier),]
      # idem 3ieme dose
      tempo_vhb_dose <- tempo_vhb_dose[duplicated(tempo_vhb_dose$nodossier),]
      tempo_vhb_3dose <- tempo_vhb_dose[!duplicated(tempo_vhb_dose$nodossier),]
      # idem 4ieme dose
      tempo_vhb_dose <- tempo_vhb_dose[duplicated(tempo_vhb_dose$nodossier),]
      tempo_vhb_4dose <- tempo_vhb_dose[!duplicated(tempo_vhb_dose$nodossier),]
      
      sortie <- list(
            ("vhb, 1iere dose, resume des ages"), summary(as.numeric(tempo_vhb_1dose$age_vacci_mois[tempo_vhb_1dose$acte_code == "P1"])),
            ("vhb, 2ieme dose, resume des ages"), summary(as.numeric(tempo_vhb_2dose$age_vacci_mois[tempo_vhb_2dose$acte_code == "P2"])),
            ("vhb, 3ieme dose, resume des ages"), summary(as.numeric(tempo_vhb_3dose$age_vacci_mois[tempo_vhb_3dose$acte_code == "P3"])),
            ("vhb, 4ieme dose, resume des ages"), summary(as.numeric(tempo_vhb_4dose$age_vacci_mois[tempo_vhb_4dose$acte_code == "R01"])),
            ("vhb, CV, 1 dose"), dim(tempo_vhb_1dose)[1] / tempo_population * 100,
            ("vhb, CV, 2 doses"), dim(tempo_vhb_2dose)[1] / tempo_population * 100,
            ("vhb, CV, 3 doses"), dim(tempo_vhb_3dose)[1] / tempo_population * 100,
            ("vhb, CV, 4 doses"), dim(tempo_vhb_4dose)[1] / tempo_population * 100
     )
      return(sortie)
}

for (i in 2010:2015) {
      print(paste("*************",i,"***************"))
      print(f_vhb(donnees,i))
}


# fin analyse hepatite b
####################################################################################################


####################################################################################################
# ANALYSE MENINGO C
############################################################################

# meningoC_2015 <- PMI_2015[PMI_2015$vaccin_code %in% c("MEN C"),]
# meningoACYW_2015 <- PMI_2015[PMI_2015$vaccin_code %in% c("MEN C","M ACYW","M ACYW135"),] 
# meningoACYW_2014 <- PMI_2014[PMI_2014$vaccin_code %in% c("MEN C","M ACYW","M ACYW135"),]
# meningoACYW_2013 <- PMI_2013[PMI_2013$vaccin_code %in% c("MEN C","M ACYW","M ACYW135"),]
# meningoACYW_2012 <- PMI_2012[PMI_2012$vaccin_code %in% c("MEN C","M ACYW","M ACYW135"),]
# meningoACYW_2011 <- PMI_2011[PMI_2011$vaccin_code %in% c("MEN C","M ACYW","M ACYW135"),]

# # nb de chaque vaccin
# table(
#       donnees$vaccin_code[
#             donnees$vaccin_code %in% c("MEN C","M ACYW","M ACYW135") 
#       ]
# )
# 
# # annee vacci meningo C
# plot(table(format(meningoACYW$dateopv,"%Y")))
# plot(table(meningoACYW$age_vacci)) # permet bien de voir une premiere coupure a 3 ans
# # et une 2ieme a 24 ans
# 
# # distribution des ages de vacci
# 
# p_meningoC <- ggplot(meningoC)
# p_meningoC + geom_bar(aes(x=age_vacci,colour = "blue"))
# 
# p_meningoACYW <- ggplot(meningoACYW)
# p_meningoACYW + geom_bar(aes(x=age_vacci,colour = "blue"))
#####################
################################################################################

# fonction sur vaccin C uniquement

f_menC <- function (df_men,an){
      ############################################################
      # preparation
      
      # restriction du df d'entree a l'annee voulue et au vaccin
      ddn_min <- as.Date(paste(an-3,"-01-01",sep="")) # on veut les enfants ayants 3/4 ans l'annee choisie
      ddn_max <- as.Date(paste(an-3,"-12-31",sep=""))
      tempo_men <- df_men[df_men$datenaiss > ddn_min &
                          df_men$datenaiss < ddn_max,]
      
      tempo_population <- comptage(tempo_men,an,3) # calcul population denominateur avant restriction sur vaccines
      tempo_men <- tempo_men[tempo_men$vaccin_code %in% menC,]
      
      # on classe suivant le n° de dossier et l'anciennete de l'opv
      tempo_men <- tempo_men[order(tempo_men$nodossier,tempo_men$dateopv),]
      
      
      #############################################################
      # debut des calculs des nombres de vaccines
      
      # creation d'un df avec une seule repetition pour chaque dossier
      tempo_men <- tempo_men[!duplicated(tempo_men$nodossier),]
      
     
      sortie <- list(
            round(dim(tempo_men)[1] / tempo_population * 100,2),
            binom.test(dim(tempo_men)[1],tempo_population)$conf.int,
            summary(as.numeric(tempo_men$age_vacci_mois))
            
     )
      return(sortie)
}

# fonction sur tout vaccin ACYW
f_menACYW <- function (df_men,an){
      
      ############################################################
      # preparation
      
      # restriction du df d'entree a l'annee voulue et au vaccin
      ddn_min <- as.Date(paste(an-3,"-01-01",sep="")) # on veut les enfants ayants 3/4 ans l'annee choisie
      ddn_max <- as.Date(paste(an-3,"-12-31",sep=""))
      tempo_men <- df_men[df_men$datenaiss > ddn_min &
                          df_men$datenaiss < ddn_max,]
      
      tempo_population <- comptage(tempo_men,an,3) # calcul population denominateur avant restriction sur vaccines
      tempo_men <- tempo_men[tempo_men$vaccin_code %in% menACYW,]
      
      # on classe suivant le n° de dossier et l'anciennete de l'opv
      tempo_men <- tempo_men[order(tempo_men$nodossier,tempo_men$dateopv),]
      
      
      #############################################################
      # debut des calculs des nombres de vaccines
      
      # creation d'un df avec une seule repetition pour chaque dossier
      tempo_men <- tempo_men[!duplicated(tempo_men$nodossier),]
      
      sortie <- list(
            round(dim(tempo_men)[1] / tempo_population * 100,2),
            binom.test(dim(tempo_men)[1],tempo_population)$conf.int,
            summary(as.numeric(tempo_men$age_vacci_mois))
     )
      return(sortie)
}
###################################################
# boucle de sortie des resultats

print("*****************************
      ****** MEN C *****************
      ******************************")
sortie <- matrix(NA,7,6)
colnames(sortie) <- 2010:2015
rownames(sortie) <- c("CV en %","IC inf (%)","IC sup (%)","1er quartile","médiane","moyenne","3ième quartile")

#remplissage du tableau de sortie
for (i in 2010:2015) {
      tempo <- f_menC(donnees,i)
      sortie[1,i-2009] <- tempo[[1]] # CV
      sortie[2,i-2009] <- tempo[[2]][1]*100 ; sortie[3,i-2009] <- tempo[[2]][2]*100 # IC
      for (j in 4:7)
            sortie[j,i-2009] <- tempo[[3]][j-2] # repartition des ages en mois
}
round(sortie,2)
rm(tempo,sortie) # on efface les fichiers temporaires pour faire de la place

print("*****************************
      ****** MEN ACYW **************
      ******************************")
sortie <- matrix(NA,7,6)
colnames(sortie) <- 2010:2015
rownames(sortie) <- c("CV en %","IC inf (%)","IC sup (%)","1er quartile","médiane","moyenne","3ième quartile")

#remplissage du tableau de sortie
for (i in 2010:2015) {
      tempo <- f_menACYW(donnees,i)
      sortie[1,i-2009] <- tempo[[1]] # CV
      sortie[2,i-2009] <- tempo[[2]][1]*100 ; sortie[3,i-2009] <- tempo[[2]][2]*100 # IC
      for (j in 4:7)
            sortie[j,i-2009] <- tempo[[3]][j-2] # repartition des ages en mois
}
round(sortie,2)
rm(tempo,sortie) # on efface les fichiers temporaires pour faire de la place
      

print("*****************************
      ****** MEN C *****************
      ******************************
      ****** PAR CANTONS ***********
      ******************************")



for (c in levels(donnees$canton)){
      print(paste("*****************************",c,"************************************************************"))
      tempo_canton <- donnees[donnees$canton == c,]
      
      sortie <- matrix(NA,7,6)
      colnames(sortie) <- 2010:2015
      rownames(sortie) <- c("CV en %","IC inf (%)","IC sup (%)","1er quartile","médiane","moyenne","3ième quartile")

      #remplissage du tableau de sortie
      for (i in 2010:2015) {
            tempo <- f_menC(tempo_canton,i)
            sortie[1,i-2009] <- tempo[[1]] # CV
            sortie[2,i-2009] <- tempo[[2]][1]*100 ; sortie[3,i-2009] <- tempo[[2]][2]*100 # IC
            for (j in 4:7)
                  sortie[j,i-2009] <- tempo[[3]][j-2] # repartition des ages en mois
      }
      print(round(sortie,2))
      rm(tempo,tempo_canton,sortie) # on efface les fichiers temporaires pour faire de la place
      
}

print("*****************************
      ****** MEN ACYW **************
      ******************************
      ****** PAR CANTONS ***********
      ******************************")



for (c in levels(donnees$canton)){
      print(paste("*****************************",c,"************************************************************"))
      tempo_canton <- donnees[donnees$canton == c,]
      
      sortie <- matrix(NA,7,6)
      colnames(sortie) <- 2010:2015
      rownames(sortie) <- c("CV en %","IC inf (%)","IC sup (%)","1er quartile","médiane","moyenne","3ième quartile")

      #remplissage du tableau de sortie
      for (i in 2010:2015) {
            tempo <- f_menACYW(tempo_canton,i)
            sortie[1,i-2009] <- tempo[[1]] # CV
            sortie[2,i-2009] <- tempo[[2]][1]*100 ; sortie[3,i-2009] <- tempo[[2]][2]*100 # IC
            for (j in 4:7)
                  sortie[j,i-2009] <- tempo[[3]][j-2] # repartition des ages en mois
      }
      print(round(sortie,2))
      rm(tempo,tempo_canton,sortie) # on efface les fichiers temporaires pour faire de la place
      
}

# fin analyse meningo C
####################################################################################################


# ####################################################################################################
# # ANALYSE ROR
# #################################################################################
# 
# ################
# ##########################################################################################
# 
# 
# # creation de la fonction calcul CV pour tempo_ror
# f_ror <- function (df_ror,an){
#       ############################################################
#       # preparation
#       
#       # restriction du df d'entree a l'annee voulue et au vaccin
#       ddn_min <- as.Date(paste(an-3,"-01-01",sep="")) # on veut les enfants ayants 3/4 ans l'annee choisie
#       ddn_max <- as.Date(paste(an-3,"-12-31",sep=""))
#       tempo_ror <- df_ror[df_ror$datenaiss > ddn_min &
#                           df_ror$datenaiss < ddn_max,]
#       
#       tempo_population <- comptage(tempo_ror,an,3) # calcul population denominateur avant restriction sur vaccines
#       tempo_ror <- tempo_ror[tempo_ror$vaccin_code %in% ror,]
#       
#       # on classe suivant le n° de dossier et l'anciennete de l'opv
#       tempo_ror <- tempo_ror[order(tempo_ror$nodossier,tempo_ror$dateopv),]
#       
#       #############################################################
#       # debut des calculs des nombres de vaccines
#       
#       # creation d'un df avec une seule repetition pour chaque dossier => nb de vacci 1 dose
#       tempo_ror_1dose <- tempo_ror[!duplicated(tempo_ror$nodossier),]
#       # pour 2 doses, on prend les autres vacci donc, sans "!", puis on réapplique
#       # utilisation d'un df "pivot" tempo_ror_dose
#       tempo_ror_dose <- tempo_ror[duplicated(tempo_ror$nodossier),] # toutes les vacci sauf la premiere dose
#       tempo_ror_2dose <- tempo_ror_dose[!duplicated(tempo_ror_dose$nodossier),]
#       # idem 3ieme dose
#       tempo_ror_dose <- tempo_ror_dose[duplicated(tempo_ror_dose$nodossier),]
#       tempo_ror_3dose <- tempo_ror_dose[!duplicated(tempo_ror_dose$nodossier),]
#       # idem 4ieme dose
#       tempo_ror_dose <- tempo_ror_dose[duplicated(tempo_ror_dose$nodossier),]
#       tempo_ror_4dose <- tempo_ror_dose[!duplicated(tempo_ror_dose$nodossier),]
#       
#       sortie1 <- list(
#             round(dim(tempo_ror_1dose)[1] / tempo_population * 100,2),
#             binom.test(dim(tempo_ror_1dose)[1],tempo_population)$conf.int,
#             summary(as.numeric(tempo_ror_1dose$age_vacci_mois[tempo_ror_1dose$acte_code == "P1"]))
#             )
#       sortie2 <- list(
#             round(dim(tempo_ror_2dose)[1] / tempo_population * 100,2),
#             binom.test(dim(tempo_ror_2dose)[1],tempo_population)$conf.int,
#             summary(as.numeric(tempo_ror_2dose$age_vacci_mois[tempo_ror_2dose$acte_code == "P2"]))
#             )
#       sortie3 <- list(
#             round(dim(tempo_ror_3dose)[1] / tempo_population * 100,2),
#             binom.test(dim(tempo_ror_3dose)[1],tempo_population)$conf.int,
#             summary(as.numeric(tempo_ror_3dose$age_vacci_mois[tempo_ror_3dose$acte_code == "P3"]))
#             )
#       sortie4 <- list(
#             round(dim(tempo_ror_4dose)[1] / tempo_population * 100,2),
#             binom.test(dim(tempo_ror_4dose)[1],tempo_population)$conf.int,
#             summary(as.numeric(tempo_ror_4dose$age_vacci_mois[tempo_ror_4dose$acte_code == "R01"]))
#             )
#       sortie <- list(sortie1,sortie2,sortie3,sortie4) # liste de listes ...
#       return(sortie)
# }
# 
# 
# ###################################################
# # boucle de sortie des resultats
# 
# print("*****************************
#       ****** ror *****************
#       ******************************")
# 
# #remplissage des tableaux de sortie
# for (d in 1:4){ # d pour dose
#       print(paste("---------------------------",d," dose(s) ----------------------------"))
#       
#       sortie <- matrix(NA,7,6)
#       colnames(sortie) <- 2010:2015
#       rownames(sortie) <- c("CV en %","IC inf (%)","IC sup (%)","1er quartile","médiane","moyenne","3ième quartile")
#    
#       for (i in 2010:2015) {
#             tempo <- f_ror(donnees,i)[[d]] # on ne recupere que la "d"ieme liste dans la sortie
#             sortie[1,i-2009] <- tempo[[1]] # CV
#             sortie[2,i-2009] <- tempo[[2]][1]*100 ; sortie[3,i-2009] <- tempo[[2]][2]*100 # IC
#             for (j in 4:7)
#                   sortie[j,i-2009] <- tempo[[3]][j-2] # repartition des ages en mois
#       }
#       print(round(sortie,2))
#       rm(tempo,sortie) # on efface les fichiers temporaires pour faire de la place
# }
# 
# 
# 
# 
# 
# print("*****************************
#       ****** ror *****************
#       ******************************
#       ****** PAR CANTONS ***********
#       ******************************")
# 
# 
# 
# for (c in levels(donnees$canton)){
#       print(paste("*****************************",c,"************************************************************"))
#       tempo_canton <- donnees[donnees$canton == c,]
#       for (d in 1:4) {
#             print(paste("---------------------------",d," dose(s) ----------------------------"))
#             sortie <- matrix(NA,7,6)
#             colnames(sortie) <- 2010:2015
#             rownames(sortie) <- c("CV en %","IC inf (%)","IC sup (%)","1er quartile","médiane","moyenne","3ième quartile")
#       
#             #remplissage du tableau de sortie
#             for (i in 2010:2015) {
#                   tempo <- f_ror(tempo_canton,i)[[d]] # on ne recupere que la "d"ieme liste dans la sortie
#                   sortie[1,i-2009] <- tempo[[1]] # CV
#                   sortie[2,i-2009] <- tempo[[2]][1]*100 ; sortie[3,i-2009] <- tempo[[2]][2]*100 # IC
#                   for (j in 4:7)
#                         sortie[j,i-2009] <- tempo[[3]][j-2] # repartition des ages en mois
#             }
#             print(round(sortie,2))
#             # rm(tempo,tempo_canton,sortie) # on efface les fichiers temporaires pour faire de la place
#       }
# }
# 
# 
# 
# 
# # interet a separer les 3 valences 
# tempo <- donnees[donnees$vaccin_code %in% ror &
#                        donnees$datenaiss > as.Date("2007-01-01") &
#                        donnees$datenaiss < as.Date("2012-12-31"),]
# 
# table(factor(tempo$vaccin_code))
# 
# # 1 rubeole seule, 2 rougeole rubeole, 8 rougeole seule et 12385 ROR, donc analyse sur ROR seul, 
# # => ror <- c("ROR")
# 
# 
# # fin analyse ror
# ####################################################################################################

















table(factor(DTPolio_2015$acte_code))
summary(DTPolio_2015$age_vacci[DTPolio_2015$acte_code == "P3"])



# verif acte_code
table(factor(coq_2015_1dose$acte_code))
DTPolio_2015_3dose$age_vacci_mois[DTPolio_2015_3dose$acte_code == "R11M"]
DTPolio_2015_2dose$vaccin_code[DTPolio_2015_2dose$acte_code != "P2"]

summary(DTPolio_2015_4dose$age_vacci_mois[DTPolio_2015_3dose$acte_code == "R01"])




#########*****************############

cantons <- read.csv2("donnees/communes_cantons.csv")
cantons$Commune <- toupper(cantons$Commune)
cantons$Commune <- factor(cantons$Commune)

sort(cantons$Commune)
summary(cantons)


difference <- sort(setdiff(donnees$commune[as.numeric(as.character(donnees$cp)) >= 20000 &
                             as.numeric(as.character(donnees$cp)) < 20200],cantons$Commune))
cantons$Commune

table(factor(donnees$commune[as.numeric(as.character(donnees$cp)) >= 20000 &
                             as.numeric(as.character(donnees$cp)) < 20200]))



extraction <- donnees[donnees$commune %in% difference,]
write.csv2(extraction,"extraction.csv",row.names = F)

table(factor(donnees$commune[donnees$cp == "20167"] ))


