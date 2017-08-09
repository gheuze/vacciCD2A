graphe <- read.csv2("sorties/dtp/dtp - pour graphe.csv")

colnames(graphe) <- 1995:2013
rownames(graphe) <- c("primovaccination à 24 mois", "rappel à 24 mois","primovaccination à 4 ans", "rappel à 4 ans")

png("sorties/dtp/evolution_CV_24-4.png",width=25 ,height=15, units="cm",res = 200)
matplot(t(graphe),type="b",pch=1:4,xlab = "années", ylim = c(75,100),ylab = "taux de couverture vaccinal en %",
         main = "Evolution de la couverture vaccinale DTP à l'âge de 24 mois et 4 ans 
         - Corse-du-Sud, 1997-2015 - axe des ordonnées débutant à 75 %", xaxt="n",
        col = couleurs)
axis(1,at=1:19,labels=1997:2015)
legend(13,80,legend = rownames(graphe), pch=1:4,col = couleurs, cex=0.8,bty="n",lwd=1)
dev.off()


reg <- as.numeric(graphe[4,-(1:10)]) ; date <- as.numeric(2005:2013)
summary(lm(reg~date))
