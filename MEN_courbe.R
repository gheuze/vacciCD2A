# comparaison niveau national, CS24 et departemental

graphe <- read.csv2("sorties/men/men_24 - pour graphe.csv",header = F)

colnames(graphe) <- 2010:2016
rownames(graphe) <- c("EGB, CV nationale","EGB, CV régionale","données Corse-du-Sud")

couleurs <- c("red","green","blue")
png("sorties/men/evolution_CV_nat-reg-PMI.png",width=25 ,height=15, units="cm",res = 400)
      matplot(t(graphe),type="b",pch = c(0,0,3),xlab = "années", ylim = c(0,100),ylab = "taux de couverture vaccinal en %",
               main = "Comparaison de l'évolution de la couverture vaccinale méningocoque C
      à l'âge de 24 mois suivant les sources, 2010-2015", xaxt="n",lty = c(1,3,3),col = couleurs)
      axis(1,at=1:7,labels=2010:2016)
      legend(5.5,20,legend = rownames(graphe), pch=c(0,0,3),col = couleurs, cex=0.8,bty="n",lwd=1,lty = c(1,3,3))
      legend(5.5,7,legend = "objectif de vaccination", pch="",col = "darkorange", cex=0.8,bty="n",lwd=1,lty = 2)
      abline(h = 95, col = "darkorange",lty = 2)
dev.off()

# comparaison 24 mois - 4 ans
graphe <- read.csv2("sorties/men/men_24-4.csv",header = T)

colnames(graphe) <- 2006:2011
rownames(graphe) <- c("CV à 24 mois", "CV à 4 ans")
couleurs <- c("blue")

png("sorties/men/evolution_CV_24-4.png",width=25 ,height=15, units="cm",res = 400)
matplot(t(graphe),type="b",pch=1:2,xlab = "année de naissance de l'enfant", ylim = c(0,100),ylab = "taux de couverture vaccinal en %",
         main = "Evolution de la couverture vaccinale méningocoque C à l'âge de 24 mois et 4 ans,
        pour les enfants nés entre 2006 et 2011", xaxt="n", col = couleurs, lty = 3)
axis(1,at=1:6,labels=2006:2011)
legend(5,17,legend = rownames(graphe), pch=1:2,col = couleurs, cex=0.8,bty="n",lwd=1,lty = 3)
legend(5,8,legend = "objectif de vaccination", pch="",col = "darkorange", cex=0.8,bty="n",lwd=1,lty = 2)
abline(h = 95, col = "darkorange",lty = 2)
dev.off()

# date <- as.numeric(1995:2013)
# for (i in 1:4) {
#       reg <- as.numeric(graphe[i,])  
#       print(summary(lm(reg~date)))
# }

