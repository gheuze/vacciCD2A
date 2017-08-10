# comparaison niveau national, CS24 et departemental

graphe <- read.csv2("sorties/vhb/vhb_24 - Copie.csv",header = T)

colnames(graphe) <- 1995:2015
rownames(graphe) <- c("vacci nat","CS24 d�p","d�pa d�part")


png("sorties/vhb/evolution_CV_nat-CS24-PMI.png",width=25 ,height=15, units="cm",res = 400)
matplot(t(graphe),type="b",pch=1:3,xlab = "ann�es", ylim = c(0,100),ylab = "taux de couverture vaccinal en %",
         main = "Comparaison de l'�volution de la couverture vaccinale coqueluche � l'�ge de 24 mois suivant les sources 
         - Corse-du-Sud, 1997-2015 - axe des ordonn�es d�butant � 75 %", xaxt="n",
        col = couleurs)
axis(1,at=1:21,labels=1995:2015)
legend(18,20,legend = rownames(graphe), pch=1:3,col = couleurs, cex=0.8,bty="n",lwd=1)
abline(h = 95, col = "darkorange",lty = 2)
dev.off()

# date <- as.numeric(1995:2013)
# for (i in 1:4) {
#       reg <- as.numeric(graphe[i,])  
#       print(summary(lm(reg~date)))
# }
# 

# comparaison 24 mois - 4 ans
graphe <- read.csv2("sorties/vhb/vhb_24-4.csv",header = F)

colnames(graphe) <- 1993:2011
rownames(graphe) <- c("vaccination � 24 mois", "vaccination � 4 ans")

png("sorties/vhb/evolution_CV_24-4.png",width=25 ,height=15, units="cm",res = 400)
matplot(t(graphe),type="b",pch=1:2,xlab = "ann�e de naissance de l'enfant", ylim = c(0,100),ylab = "taux de couverture vaccinal en %",
         main = "Evolution de la couverture vaccinale VHB � l'�ge de 24 mois et 4 ans 
         - Corse-du-Sud, pour les enfants n�s entre 1993 et 2011", xaxt="n",
        col = couleurs)
axis(1,at=1:19,labels=1993:2011)
legend(12,15,legend = rownames(graphe), pch=1:2,col = couleurs, cex=0.8,bty="n",lwd=1)
abline(h = 95, col = "darkorange",lty = 2)
dev.off()

# date <- as.numeric(1995:2013)
# for (i in 1:4) {
#       reg <- as.numeric(graphe[i,])  
#       print(summary(lm(reg~date)))
# }
