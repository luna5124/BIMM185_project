library(ggplot2)
library(plyr)

operon = read.table("operons_rho.out")

png("operon_plot.png")

ggplot(operon, aes(x = V1, y = V2)) + geom_point()+
  xlab("distance") +
  ylab("Correlation coefficient") +
  ggtitle("Correlation Coefficient VS Distance between Operon gene pairs")

dev.off()

png("operon_plot_bins.png")

operon_cut = ddply(operon, .(cut(operon$V1, 10)), colwise(mean))
ggplot(operon_cut, aes(x = operon_cut[,1], y = V2)) + geom_point() + 
  xlab("distance") +
  ylab("Correlation coefficient") +
  ggtitle("Correlation Coefficient VS Distance between Operon gene pairs")+
  theme(axis.text.x = element_blank())
dev.off()


rhos = read.table("dis_rho_new.out")

png("dis_plot.png")
ggplot(rhos, aes(x = V1, y = V2)) + geom_point()+
xlab("distance") +
ylab("Correlation coefficient") +
ggtitle("Correlation Coefficient VS Distance between TF shared gene pairs")

dev.off()

png("dis_plot_bins.png")
rhos_cut = ddply(rhos, .(cut(rhos$V1, 10)), colwise(mean))
ggplot(rhos_cut, aes(x = rhos_cut[,1], y = V2)) + geom_point() + 
xlab("distance") +
ylab("Correlation coefficient") +
ggtitle("Correlation Coefficient VS Distance between TF shared gene pairs")+
theme(axis.text.x = element_blank())

dev.off()


l = rep(4641652, nrow(rhos))

rhos = cbind(rhos, pmin(l - rhos[,1], rhos[,1]))

png("dis_plot_new.png")
ggplot(rhos, aes(x = V1, y = V2)) + geom_point()+
xlab("distance") +
ylab("Correlation coefficient") +
ggtitle("Correlation Coefficient VS Distance between TF shared gene pairs")

dev.off()


rhos_cut = ddply(rhos, .(cut(rhos[,3], 10)), colwise(mean))
png("dis_plot_bins_new.png")
ggplot(rhos_cut, aes(x = rhos_cut[,1], y = rhos_cut[,2])) + geom_point() + 
xlab("distance") +
ylab("Correlation coefficient") +
ggtitle("Correlation Coefficient VS Distance between TF shared gene pairs")+
theme(axis.text.x = element_blank())
dev.off()
