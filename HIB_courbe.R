# comparaison niveau national, CS24 et departemental

graphe <- read.csv2("sorties/hib/hib_24 - Copie.csv",header = T)

colnames(graphe) <- 1995:2015
rownames(graphe) <- c("CS24 nat","CS24 Corse-du-Sud","base départementale")

couleurs <- c("red","blue","blue")
png("sorties/hib/evolution_CV_nat-CS24-PMI.png",width=25 ,height=15, units="cm",res = 400)
matplot(t(graphe),type="b",pch=c(1,1,3),xlab = "années", ylim = c(60,100),ylab = "taux de couverture vaccinal en %",
         main = "Comparaison de l'évolution de la couverture vaccinale Hib à l'âge de 24 mois suivant les sources 
         - Corse-du-Sud, 1995-2015, début de l'axe à 60 %", xaxt="n",col = couleurs, lty =c(1,3,3))
axis(1,at=1:21,labels=1995:2015)
legend(15,70,legend = rownames(graphe), pch=c(1,1,3),col = couleurs, cex=0.8,bty="n",lwd=1,lty=c(1,3,3))
legend(15,64,legend = "objectif de vaccination", pch="",col = "darkorange", cex=0.8,bty="n",lwd=1,lty = 2)
abline(h = 95, col = "darkorange",lty = 2)
dev.off()

# date <- as.numeric(1995:2013)
# for (i in 1:4) {
#       reg <- as.numeric(graphe[i,])  
#       print(summary(lm(reg~date)))
# }
# 

# comparaison 24 mois - 4 ans
graphe <- read.csv2("sorties/hib/hib_24-4.csv",header = T)

colnames(graphe) <- 1993:2011
rownames(graphe) <- c("vaccination à 24 mois", "vaccination à 4 ans")

couleurs <- c("blue")
png("sorties/hib/evolution_CV_24-4.png",width=25 ,height=15, units="cm",res = 400)
matplot(t(graphe),type="b",pch=3:2,xlab = "année de naissance de l'enfant", ylim = c(60,100),ylab = "taux de couverture vaccinal en %",
         main = "Evolution de la couverture vaccinale Hib à l'âge de 24 mois et 4 ans 
         - Corse-du-Sud, pour les enfants nés entre 1993 et 2011, début de l'axe à 60 %", xaxt="n",
        col = couleurs, lty =c(1,3))
axis(1,at=1:19,labels=1993:2011)
legend(15,68,legend = rownames(graphe), pch=3:2,col = couleurs, cex=0.8,bty="n",lwd=1, lty =c(1,3))
legend(15,64,legend = "objectif de vaccination", pch="",col = "darkorange", cex=0.8,bty="n",lwd=1,lty = 2)
abline(h = 95, col = "darkorange",lty = 2)
dev.off()

# date <- as.numeric(1995:2013)
# for (i in 1:4) {
#       reg <- as.numeric(graphe[i,])  
#       print(summary(lm(reg~date)))
# }

