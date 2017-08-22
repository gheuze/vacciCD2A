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


