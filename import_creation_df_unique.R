##################################
#
#   analyse donnees vacci CD 2A -- import et creation d'un df unique
#
# V1  août 2017
##################################


rm(list=ls())
setwd("I:\\DIRECTION ACTION TERRITORIALE DE SANTE\\VSSE\\CIRE\\Dossiers MI\\Vaccination\\etude_CG2A\\analyse")

# import donnees administratives
donnees1 <- read.csv2("donnees/export_1.csv")

donnees1$nodossier <- factor(donnees1$nodossier)
donnees1$datenaiss <- as.Date(donnees1$datenaiss, format="%d/%m/%Y")

# import donnees medicales
donnees2 <- read.csv2("donnees/export_2.csv")
donnees2$noetcv <- factor(donnees2$noetcv)
donnees2$dateopv <- as.Date(donnees2$dateopv,format="%d/%m/%Y")
donnees2$datemin <- as.Date(donnees2$datemin,format="%d/%m/%Y")
donnees2$datetyp <- as.Date(donnees2$datetyp,format="%d/%m/%Y")
donnees2$datemax <- as.Date(donnees2$datemax,format="%d/%m/%Y")



# fusion des 2 df et supression des df d'origine
donnees <- merge(donnees1,donnees2,by.x="nodossier",by.y="noetcv")
rm(donnees1,donnees2)

# analyse des donnees jusqu'au 31/12/2015 car base 2016 non complete
donnees <- donnees[donnees$dateopv <= as.Date("2015-12-31"),]

# liste des vaccins ne nous interesssant pas, on enleve pour alleger la base
a_exclure <- c("BCG","ENC JAP","FJ","G","L","RA","SHA","SHB","TAB","HPV",
               "TVI","TVI + HA","V","ZONA")

donnees <- donnees[!donnees$vaccin_code %in% a_exclure,]
donnees$vaccin_code <- factor(donnees$vaccin_code)

# calcul d'un age a la vacci
donnees$age_vacci <- as.numeric(floor((donnees$dateopv-donnees$datenaiss)/365.25))

# calcul d'un age en mois
donnees$age_vacci_mois <- difftime(donnees$dateopv, donnees$datenaiss, units = "days")
donnees$age_vacci_mois <- donnees$age_vacci_mois/365.25
donnees$age_vacci_mois <- round(donnees$age_vacci_mois*12,0)

# on a mis toutes les dates de naissance au 15 du mois pour anonymisation, 
# donc age_vacci vaut 0 pour ceux qui ont un age == -1. Par contre, une personne avec age_vacci à -39, 
# donc on enleve cette ligne
donnees <- donnees[donnees$age_vacci >= -1,] # on enlève celui qui a un age inferieur a -1
donnees$age_vacci[donnees$age_vacci == -1] <- 0 # on met a 0, tout ceux qui ont un age a -1

# attribution des cantons par CP, si celui-ci unique
donnees$canton[donnees$cp %in% c("20123","20128","20132","20134","20138","20140","20142","20148",
                                 "20153","20157","20166","20168","20173","20190")] <- "Taravo-Ornano"
donnees$canton[donnees$cp %in% c("20111","20115","20121","20125","20126","20130","20139","20141",
                                 "20147","20150","20151","20160")] <- "Sevi-Sorru-Cinarca"
donnees$canton[donnees$cp %in% c("20100","20110","20112","20113","20116","20122","20127","20143",
                                 "20152","20164","20165")] <- "Sartenais Valinco"
donnees$canton[donnees$cp %in% c("20114","20124","20131","20135","20137","20145","20146","20169",
                                 "20170","20171")] <- "Bavella-Grand sud"
donnees$canton[donnees$cp %in% c("20119","20133","20136","20163","20172")] <- "Gravona-Prunelli"
donnees$canton[donnees$cp %in% c("20000","20129")] <- "Ajaccio"

# passage par les communes pour les CP sur plusiseurs cantons
donnees$canton[donnees$commune %in% c("CAURO","SAINT JEAN DE PISCIATELLO","ECCICA SUARELLA")] <- 
      "Taravo-Ornano"
donnees$canton[donnees$commune %in% c("TOLLA","OCANA","AFA","CORTICCHIATO","PISCIA ROSSA","CUTTOLI",
                                      "SARROLA","VOLPAJA","CUTTOLI CORTICCHIATO","SARROLA CARCOPINO",
                                      "APPIETTO","TAVACO","PERI","VALLE DI MEZZANA")] <- 
      "Gravona-Prunelli"
donnees$canton[donnees$commune %in% c("VILLANOVA","AJACCIO","ALATA","MEZZAVIA","CARDIGLIONE")] <- 
      "Ajaccio"
donnees$canton <- factor(donnees$canton)

save.image("donnees/df_pour_analyse.RData")


# fin creation d'un df unique
##########################################################