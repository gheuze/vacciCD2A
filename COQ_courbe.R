
# comparaison niveau national, CS24 et departecoqtal

graphe <- read.csv2("sorties/coq/coq - pour graphe.csv",header = T)

colnames(graphe) <- 1995:2015
rownames(graphe) <- c("CV nationale issue des CS 24", "CV nationale, issue des CS 24",
                      "CV 2A issue des CS 24", "CV 2A, issue des CS 24",
                      "CV 2A issue de le base", "CV 2A, issue de le base")

couleurs <- c("red",rep("blue",2))

png("sorties/coq/evolution_CV_nat-dep-PMI%02d.png",width=25 ,height=15, units="cm",res = 400)
# figure sur les primovacci, cad lignes impaires du df, soit graphe[c(T,F),]
      matplot(t(graphe[c(T,F),]),type = "b",pch = c(rep(1,2),3), xlab = "années", ylim = c(75,100),ylab = "taux de couverture vaccinal en %",
               main = "Comparaison de l'évolution de la primovaccination coqueluche
      à l'âge de 24 mois suivant les sources, 1995-2015 - axe des ordonnées débutant à 75 %", xaxt="n",lty = c(1,rep(3,2)),
              col = couleurs)
      axis(1,at=1:21,labels=1995:2015)
      legend(15,80,legend = rownames(graphe[c(T,F),]), pch = c(1,1,3), col = couleurs, cex=0.8,bty="n",lwd=1,lty = c(1,rep(3,2)))
      legend(15,77,legend = "objectif de vaccination", pch="",col = "darkorange", cex=0.8,bty="n",lwd=1,lty = 2)
      abline(h = 95, col = "darkorange",lty = 2)

# figure sur les rappels, cad lignes paires du df, soit graphe[c(F,T),]
      matplot(t(graphe[c(F,T),]),type = "b",pch = c(rep(5,2),4), xlab = "années", ylim = c(75,100),ylab = "taux de couverture vaccinal en %",
               main = "Comparaison de l'évolution du rappel coqueluche
      à l'âge de 24 mois suivant les sources, 1995-2015 - axe des ordonnées débutant à 75 %", xaxt="n",lty = c(1,rep(3,2)),
              col = couleurs)
      axis(1,at=1:21,labels=1995:2015)
      legend(15,80,legend = rownames(graphe[c(F,T),]), pch = c(5,5,4), col = couleurs, cex=0.8,bty="n",lwd=1,lty = c(1,rep(3,2)))
      legend(15,77,legend = "objectif de vaccination", pch="",col = "darkorange", cex=0.8,bty="n",lwd=1,lty = 2)
      abline(h = 95, col = "darkorange",lty = 2)
dev.off()

#########################################################################################################################

# comparaison 24 mois - 4 ans
graphe <- read.csv2("sorties/coq/coq_24-4.csv",header = T)

colnames(graphe) <- 1993:2011
rownames(graphe) <- c("primovaccination à 24 mois", "rappel à 24 mois","primovaccination à 4 ans", "rappel à 4 ans")
couleurs <- c("blue")

png("sorties/coq/evolution_CV_24-4.png",width=25 ,height=15, units="cm",res = 400)
      matplot(t(graphe),type="b",pch = c(3,4,2,6), xlab = "année de naissance de l'enfant", ylim = c(45,100),ylab = "taux de couverture vaccinal en %",
               main = "Evolution de la couverture vaccinale coqueluche à l'âge de 24 mois et 4 ans,
              pour les enfants nés entre 1993 et 2011 - axe des ordonnées débutant à 45 %", xaxt="n", col = couleurs, lty = 3)
      axis(1,at=1:19,labels=1993:2011)
      legend(15,57,legend = rownames(graphe), pch = c(3,4,2,6),col = couleurs, cex=0.8,bty="n",lwd=1,lty = 3)
      legend(15,48,legend = "objectif de vaccination", pch="",col = "darkorange", cex=0.8,bty="n",lwd=1,lty = 2)
      abline(h = 95, col = "darkorange",lty = 2)
dev.off()