# comparaison niveau national, CS24 et departemental

graphe <- read.csv2("sorties/vhb/vhb - pour graphe.csv",header = T)

colnames(graphe) <- 1995:2015
rownames(graphe) <- c("CV nationale, issue des CS 24", "CV 2A, issue des CS 24", "CV 2A, issue de le base")

couleurs <- c("red",rep("blue",2))

png("sorties/vhb/evolution_CV_nat-dep-PMI%02d.png",width=25 ,height=15, units="cm",res = 400)
      matplot(t(graphe),type = "b",pch = c(rep(1,2),3), xlab = "années", ylim = c(0,100),ylab = "taux de couverture vaccinal en %",
               main = "Comparaison de l'évolution de la vaccination 3 doses VHB
      à l'âge de 24 mois suivant les sources, 1995-2015", xaxt="n",lty = c(1,rep(3,2)),
              col = couleurs)
      axis(1,at=1:21,labels=1995:2015)
      legend(16,25,legend = rownames(graphe), pch = c(1,1,3), col = couleurs, cex=0.8,bty="n",lwd=1,lty = c(1,rep(3,2)))
      legend(16,10.5,legend = "objectif de vaccination", pch="",col = "darkorange", cex=0.8,bty="n",lwd=1,lty = 2)
      abline(h = 95, col = "darkorange",lty = 2)
dev.off()



# comparaison 24 mois - 4 ans
graphe <- read.csv2("sorties/vhb/vhb_24-4.csv",header = T)

colnames(graphe) <- 1993:2011
rownames(graphe) <- c("vaccination à 24 mois", "vaccination à 4 ans")

couleurs <- c("blue")

png("sorties/vhb/evolution_CV_24-4.png",width=25 ,height=15, units="cm",res = 400)
      matplot(t(graphe),type = "b",pch = c(3,2),xlab = "année de naissance de l'enfant", ylim = c(0,100), 
              ylab = "taux de couverture vaccinal en %",
               main = "Evolution de la couverture vaccinale VHB à l'âge de 24 mois et 4 ans 
               - Corse-du-Sud, pour les enfants nés entre 1993 et 2011", xaxt="n",
              col = couleurs, lty = 3)
      axis(1,at=1:19,labels=1993:2011)
      legend(15,20,legend = rownames(graphe), pch = c(3,2),col = couleurs, cex=0.8,bty="n",lwd=1,lty = 3)
      legend(15,11,legend = "objectif de vaccination", pch="",col = "darkorange", cex=0.8,bty="n",lwd=1,lty = 2)
      abline(h = 95, col = "darkorange",lty = 2)
dev.off()
