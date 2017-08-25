# comparaison niveau national, CS24 et departepnetal

graphe <- read.csv2("sorties/pne/pne_24 - pour graphe.csv",header = F)

colnames(graphe) <- 2006:2015
rownames(graphe) <- c("CS24 nationale","CS24 Corse-du-Sud","données Corse-du-Sud")

couleurs <- c("red","blue","blue")
png("sorties/pne/evolution_CV_nat-reg-PMI.png",width=25 ,height=15, units="cm",res = 400)
      matplot(t(graphe),type = "b",pch = c(1,1,3),xlab = "années", ylim = c(0,100), ylab = "taux de couverture vaccinal en %",
               main = "Comparaison de l'évolution de la couverture vaccinale pneumocoque
      à l'âge de 24 mois suivant les sources, 2006-2015", xaxt="n",lty = c(1,3,3),col = couleurs)
      axis(1,at=1:10,labels=2006:2015)
      legend(8,20,legend = rownames(graphe), pch = c(1,1,3),col = couleurs, cex=0.8,bty="n",lwd=1,lty = c(1,3,3))
      legend(8,7,legend = "objectif de vaccination", pch="",col = "darkorange", cex=0.8,bty="n",lwd=1,lty = 2)
      abline(h = 95, col = "darkorange",lty = 2)
dev.off()

# comparaison 24 mois - 4 ans
graphe <- read.csv2("sorties/pne/pne_24-4.csv",header = T)

colnames(graphe) <- 2003:2011
rownames(graphe) <- c("CV à 24 mois", "CV à 4 ans")
couleurs <- c("blue")

png("sorties/pne/evolution_CV_24-4.png",width=25 ,height=15, units="cm",res = 400)
      matplot(t(graphe),type = "b",pch = c(3,2),xlab = "année de naissance de l'enfant", ylim = c(0,100),ylab = "taux de couverture vaccinal en %",
               main = "Evolution de la couverture vaccinale pneumocoque à l'âge de 24 mois et 4 ans,
              pour les enfants nés entre 2003 et 2011", xaxt="n", col = couleurs, lty = 3)
      axis(1,at=1:9,labels=2003:2011)
      legend(7,17,legend = rownames(graphe), pch = c(3,2),col = couleurs, cex = 0.8,bty = "n",lwd = 1,lty = 3)
      legend(7,8,legend = "objectif de vaccination", pch = "",col = "darkorange", cex=0.8,bty="n",lwd=1,lty = 2)
      abline(h = 95, col = "darkorange",lty = 2)
dev.off()


