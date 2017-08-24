# comparaison niveau national, CS24 et departerortal

graphe <- read.csv2("sorties/ror/ror_24 - pour graphe.csv",header = T)

colnames(graphe) <- 1995:2015
rownames(graphe) <- c("CV nationale 1 dose, issue des CS 24", "CV nationale 2 doses, issue des CS 24",
                      "CV départementale 1 dose, issue des CS 24", "CV départementale 2 doses, issue des CS 24",
                      "CV départementale 1 dose, issue de la base", "CV départementale 2 doses, issue de la base")

couleurs <- c("red",rep("blue",2))

png("sorties/ror/evolution_CV_nat-dep-PMI%02d.png",width=25 ,height=15, units="cm",res = 400)
# figure sur les 1 dose, cad lignes impaires du df, soit graphe[c(T,F),]
      matplot(t(graphe[c(T,F),]),type = "b",pch = c(rep(1,2),3), xlab = "années", ylim = c(60,100),ylab = "taux de couverture vaccinal en %",
               main = "Comparaison de l'évolution de la vaccination 1dose ROR
      à l'âge de 24 mois suivant les sources, 1995-2015 - axe des ordonnées débutant à 60 %", xaxt="n",lty = c(1,rep(3,2)),
              col = couleurs)
      axis(1,at=1:21,labels=1995:2015)
      legend(13,67,legend = rownames(graphe[c(T,F),]), pch = c(1,1,3), col = couleurs, cex=0.8,bty="n",lwd=1,lty = c(1,rep(3,2)))
      legend(13,62,legend = "objectif de vaccination", pch="",col = "darkorange", cex=0.8,bty="n",lwd=1,lty = 2)
      abline(h = 95, col = "darkorange",lty = 2)

# figure sur les 2 doses, cad lignes paires du df, soit graphe[c(F,T),]
      matplot(t(graphe[c(F,T),]),type = "b",pch = c(rep(5,2),4), xlab = "années", ylim = c(0,100),ylab = "taux de couverture vaccinal en %",
               main = "Comparaison de l'évolution de la vaccination 2 doses ROR
      à l'âge de 24 mois suivant les sources, 1995-2015", xaxt="n",lty = c(1,rep(3,2)),
              col = couleurs)
      axis(1,at=1:21,labels=1995:2015)
      legend(13,15,legend = rownames(graphe[c(F,T),]), pch = c(5,5,4), col = couleurs, cex=0.8,bty="n",lwd=1,lty = c(1,rep(3,2)))
      legend(13,3,legend = "objectif de vaccination", pch="",col = "darkorange", cex=0.8,bty="n",lwd=1,lty = 2)
      abline(h = 95, col = "darkorange",lty = 2)
dev.off()

#########################################################################################################################
# comparaison 24 mois - 4 ans
graphe <- read.csv2("sorties/ror/ror_24-4.csv",header = T)

colnames(graphe) <- 1993:2011
rownames(graphe) <- c("vaccination 1 dose à 24 mois", "vaccination 2 doses à 24 mois",
                      "vaccination 1 dose à 4 ans", "vaccination 2 doses à 4 ans")
couleurs <- c("blue")

png("sorties/ror/evolution_CV_24-4.png",width=25 ,height=15, units="cm",res = 400)
      matplot(t(graphe),type="b",pch = c(3,4,2,6), xlab = "année de naissance de l'enfant", ylim = c(0,100),ylab = "taux de couverture vaccinal en %",
               main = "Evolution de la couverture vaccinale ROR à l'âge de 24 mois et 4 ans,
              pour les enfants nés entre 1993 et 2011", xaxt="n", col = couleurs, lty = 3)
      axis(1,at=1:19,labels=1993:2011)
      legend(14,18.5,legend = rownames(graphe), pch = c(3,4,2,6),col = couleurs, cex=0.8,bty="n",lwd=1,lty = 3)
      legend(14,3,legend = "objectif de vaccination", pch="",col = "darkorange", cex=0.8,bty="n",lwd=1,lty = 2)
      abline(h = 95, col = "darkorange",lty = 2)
dev.off()
