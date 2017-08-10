graphe <- read.csv2("sorties/coq/coq - pour graphe.csv",header = F)

colnames(graphe) <- 1995:2013
rownames(graphe) <- c("primovaccination à 24 mois", "rappel à 24 mois","primovaccination à 4 ans", "rappel à 4 ans")

png("sorties/coq/evolution_CV_24-4.png",width=25 ,height=15, units="cm",res = 200)
matplot(t(graphe),type="b",pch=1:4,xlab = "années", ylim = c(75,100),ylab = "taux de couverture vaccinal en %",
         main = "Evolution de la couverture vaccinale coqueluche à l'âge de 24 mois et 4 ans 
         - Corse-du-Sud, 1997-2015 - axe des ordonnées débutant à 75 %", xaxt="n",
        col = couleurs)
axis(1,at=1:19,labels=1997:2015)
legend(12,80,legend = rownames(graphe), pch=1:4,col = couleurs, cex=0.8,bty="n",lwd=1)
abline(h = 95, col = "darkorange",lty = 2)
dev.off()

date <- as.numeric(1995:2013)
for (i in 1:4) {
      reg <- as.numeric(graphe[i,])  
      print(summary(lm(reg~date)))
}

