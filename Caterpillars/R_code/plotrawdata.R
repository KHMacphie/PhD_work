#load Dataframe.R script rm(list=ls()) setwd('/Users/s1205615/')

cater14 <- subset(cater_habitat, year=="2014", select = c(site, date, caterpillars))
cater15 <- subset(cater_habitat, year=="2015", select = c(site, date, caterpillars))
cater16 <- subset(cater_habitat, year=="2016", select = c(site, date, caterpillars))
cater17 <- subset(cater_habitat, year=="2017", select = c(site, date, caterpillars))
cater18 <- subset(cater_habitat, year=="2018", select = c(site, date, caterpillars))
cater19 <- subset(cater_habitat, year=="2019", select = c(site, date, caterpillars))
cater20 <- subset(cater_habitat, year=="2020", select = c(site, date, caterpillars))

means <- data.frame(site=cater_habitat$site, date=cater_habitat$date, cater=cater_habitat$caterpillars)
means$siteday <- paste(means$site, means$date)
means$date <- as.factor(means$date)
sitemeans <- aggregate(means$cater~means$siteday, FUN=mean)
sitemeans$siteday <- sitemeans$`means$siteday`
sitemeans$cater <- sitemeans$`means$cater`
store <- pmatch(sitemeans$siteday, means$siteday)
sitemeans <- cbind(sitemeans, means$date[store])
sitemeans$date <- sitemeans$`means$date[store]`
sitemeans <- aggregate(sitemeans$cater~sitemeans$date, FUN=mean)

means20 <- data.frame(site=cater20$site, date=cater20$date, cater=cater20$caterpillars)
means20$siteday <- paste(means20$site, means20$date)
means20$date <- as.factor(means20$date)
sitemeans20 <- aggregate(means20$cater~means20$siteday, FUN=mean)
sitemeans20$siteday <- sitemeans20$`means20$siteday`
sitemeans20$cater <- sitemeans20$`means20$cater`
store <- pmatch(sitemeans20$siteday, means20$siteday)
sitemeans20 <- cbind(sitemeans20, means20$date[store])
sitemeans20$date <- sitemeans20$`means20$date[store]`
sitemeans20 <- aggregate(sitemeans20$cater~sitemeans20$date, FUN=mean)

means19 <- data.frame(site=cater19$site, date=cater19$date, cater=cater19$caterpillars)
means19$siteday <- paste(means19$site, means19$date)
means19$date <- as.factor(means19$date)
sitemeans19 <- aggregate(means19$cater~means19$siteday, FUN=mean)
sitemeans19$siteday <- sitemeans19$`means19$siteday`
sitemeans19$cater <- sitemeans19$`means19$cater`
store <- pmatch(sitemeans19$siteday, means19$siteday)
sitemeans19 <- cbind(sitemeans19, means19$date[store])
sitemeans19$date <- sitemeans19$`means19$date[store]`
sitemeans19 <- aggregate(sitemeans19$cater~sitemeans19$date, FUN=mean)

means18 <- data.frame(site=cater18$site, date=cater18$date, cater=cater18$caterpillars)
means18$siteday <- paste(means18$site, means18$date)
means18$date <- as.factor(means18$date)
sitemeans18 <- aggregate(means18$cater~means18$siteday, FUN=mean)
sitemeans18$siteday <- sitemeans18$`means18$siteday`
sitemeans18$cater <- sitemeans18$`means18$cater`
store <- pmatch(sitemeans18$siteday, means18$siteday)
sitemeans18 <- cbind(sitemeans18, means18$date[store])
sitemeans18$date <- sitemeans18$`means18$date[store]`
sitemeans18 <- aggregate(sitemeans18$cater~sitemeans18$date, FUN=mean)

means17 <- data.frame(site=cater17$site, date=cater17$date, cater=cater17$caterpillars)
means17$siteday <- paste(means17$site, means17$date)
means17$date <- as.factor(means17$date)
sitemeans17 <- aggregate(means17$cater~means17$siteday, FUN=mean)
sitemeans17$siteday <- sitemeans17$`means17$siteday`
sitemeans17$cater <- sitemeans17$`means17$cater`
store <- pmatch(sitemeans17$siteday, means17$siteday)
sitemeans17 <- cbind(sitemeans17, means17$date[store])
sitemeans17$date <- sitemeans17$`means17$date[store]`
sitemeans17 <- aggregate(sitemeans17$cater~sitemeans17$date, FUN=mean)

means16 <- data.frame(site=cater16$site, date=cater16$date, cater=cater16$caterpillars)
means16$siteday <- paste(means16$site, means16$date)
means16$date <- as.factor(means16$date)
sitemeans16 <- aggregate(means16$cater~means16$siteday, FUN=mean)
sitemeans16$siteday <- sitemeans16$`means16$siteday`
sitemeans16$cater <- sitemeans16$`means16$cater`
store <- pmatch(sitemeans16$siteday, means16$siteday)
sitemeans16 <- cbind(sitemeans16, means16$date[store])
sitemeans16$date <- sitemeans16$`means16$date[store]`
sitemeans16 <- aggregate(sitemeans16$cater~sitemeans16$date, FUN=mean)

means15 <- data.frame(site=cater15$site, date=cater15$date, cater=cater15$caterpillars)
means15$siteday <- paste(means15$site, means15$date)
means15$date <- as.factor(means15$date)
sitemeans15 <- aggregate(means15$cater~means15$siteday, FUN=mean)
sitemeans15$siteday <- sitemeans15$`means15$siteday`
sitemeans15$cater <- sitemeans15$`means15$cater`
store <- pmatch(sitemeans15$siteday, means15$siteday)
sitemeans15 <- cbind(sitemeans15, means15$date[store])
sitemeans15$date <- sitemeans15$`means15$date[store]`
sitemeans15 <- aggregate(sitemeans15$cater~sitemeans15$date, FUN=mean)

means14 <- data.frame(site=cater14$site, date=cater14$date, cater=cater14$caterpillars)
means14$siteday <- paste(means14$site, means14$date)
means14$date <- as.factor(means14$date)
sitemeans14 <- aggregate(means14$cater~means14$siteday, FUN=mean)
sitemeans14$siteday <- sitemeans14$`means14$siteday`
sitemeans14$cater <- sitemeans14$`means14$cater`
store <- pmatch(sitemeans14$siteday, means14$siteday)
sitemeans14 <- cbind(sitemeans14, means14$date[store])
sitemeans14$date <- sitemeans14$`means14$date[store]`
sitemeans14 <- aggregate(sitemeans14$cater~sitemeans14$date, FUN=mean)

par(mfrow=c(1,1))
plot(as.numeric(as.character(sitemeans19$`sitemeans19$date`)),sitemeans19$`sitemeans19$cater`, pch=20, xlab="Ordinal Date", ylab="Mean caterpillar abundance", xlim=c(117,176))
points(as.numeric(as.character(sitemeans18$`sitemeans18$date`)),sitemeans18$`sitemeans18$cater`, pch=20, col=2)
points(as.numeric(as.character(sitemeans17$`sitemeans17$date`)),sitemeans17$`sitemeans17$cater`, pch=20, col=3)
points(as.numeric(as.character(sitemeans16$`sitemeans16$date`)),sitemeans16$`sitemeans16$cater`, pch=20, col=4)
points(as.numeric(as.character(sitemeans15$`sitemeans15$date`)),sitemeans15$`sitemeans15$cater`, pch=20, col=5)
points(as.numeric(as.character(sitemeans14$`sitemeans14$date`)),sitemeans14$`sitemeans14$cater`, pch=20, col=6)
points(as.numeric(as.character(sitemeans$`sitemeans$date`)),sitemeans$`sitemeans$cater`, pch=20, col=7)


par(mfrow=c(3,2),mar=c(4.1, 4.0, 2.5, 1.5), cex=0.9, cex.lab=1.25, cex.axis=1.2)
plot(as.numeric(as.character(sitemeans14$`sitemeans14$date`)),sitemeans14$`sitemeans14$cater`, pch=20, xlab="", ylab="", xlim=c(117,176), main="2014", ylim=c(0,0.65), yaxt='n')
axis(side = 2, at = c(0,.15, .3,.45, .6), labels = c( '0.0','', '0.3','', '0.6'))
legend("topleft", legend="n = 1106", bty = "n")
plot(as.numeric(as.character(sitemeans15$`sitemeans15$date`)),sitemeans15$`sitemeans15$cater`, pch=20, xlab="", ylab="", xlim=c(117,176), main="2015", ylim=c(0,0.22), yaxt='n')
axis(side = 2, at = c(0,.05, .1,.15, .2), labels = c( '0.0','', '0.1','', '0.2'))
legend("topleft", legend="n = 2632", bty = "n")
plot(as.numeric(as.character(sitemeans16$`sitemeans16$date`)),sitemeans16$`sitemeans16$cater`, pch=20, xlab="", ylab="Mean caterpillar abundance", xlim=c(117,176), main="2016", ylim=c(0,0.42), yaxt='n')
axis(side = 2, at = c(0,.1, .2,.3, .4), labels = c( '0.0','', '0.2','', '0.4'))
legend("topleft", legend="n = 2351", bty = "n")
plot(as.numeric(as.character(sitemeans17$`sitemeans17$date`)),sitemeans17$`sitemeans17$cater`, pch=20, xlab="", ylab="", xlim=c(117,176), main="2017", ylim=c(0,0.23), yaxt='n')
axis(side = 2, at = c(0,.05, .1,.15, .2), labels = c( '0.0','', '0.1','', '0.2'))
legend("topleft", legend="n = 6877", bty = "n")
plot(as.numeric(as.character(sitemeans18$`sitemeans18$date`)),sitemeans18$`sitemeans18$cater`, pch=20, xlab="Ordinal Date", ylab="", xlim=c(117,176), main="2018", ylim=c(0,0.82), yaxt='n')
axis(side = 2, at = c(0,.2, .4,.6, .8), labels = c( '0.0','', '0.4','', '0.8'))
legend("topleft", legend="n = 6791", bty = "n")
plot(as.numeric(as.character(sitemeans19$`sitemeans19$date`)),sitemeans19$`sitemeans19$cater`, pch=20, xlab="Ordinal Date", ylab="", xlim=c(117,176), main="2019", ylim=c(0,3), yaxt='n')
axis(side = 2, at = c(0,0.75,1.5,2.25, 3), labels = c( '0.0','', '1.5','', '3.0'))
legend("topleft", legend="n = 8310", bty = "n") # saved as 8"x8"

par(mfrow=c(4,2),mar=c(4.1, 4.0, 2.5, 1.5), cex=0.9, cex.lab=1.25, cex.axis=1.2)
plot(as.numeric(as.character(sitemeans14$`sitemeans14$date`)),sitemeans14$`sitemeans14$cater`, pch=20, xlab="", ylab="", xlim=c(117,176), main="2014", ylim=c(0,2.7), xaxt="n")
text(125,2.4, "n=1106")
axis(side = 1, at = c(120,130,140,150,160,170), labels = c( '120','','','','','170'))
plot(as.numeric(as.character(sitemeans15$`sitemeans15$date`)),sitemeans15$`sitemeans15$cater`, pch=20, xlab="", ylab="", xlim=c(117,176), main="2015", ylim=c(0,2.7), xaxt="n")
text(125,2.4, "n = 2632")
axis(side = 1, at = c(120,130,140,150,160,170), labels = c( '120','','','','','170'))
plot(as.numeric(as.character(sitemeans16$`sitemeans16$date`)),sitemeans16$`sitemeans16$cater`, pch=20, xlab="", ylab="Mean caterpillar abundance", xlim=c(117,176), main="2016", ylim=c(0,2.7), xaxt="n")
text(125,2.4, "n = 2351")
axis(side = 1, at = c(120,130,140,150,160,170), labels = c( '120','','','','','170'))
plot(as.numeric(as.character(sitemeans17$`sitemeans17$date`)),sitemeans17$`sitemeans17$cater`, pch=20, xlab="", ylab="", xlim=c(117,176), main="2017", ylim=c(0,2.7), xaxt="n")
text(125,2.4, "n = 6877")
axis(side = 1, at = c(120,130,140,150,160,170), labels = c( '120','','','','','170'))
plot(as.numeric(as.character(sitemeans18$`sitemeans18$date`)),sitemeans18$`sitemeans18$cater`, pch=20, xlab="", ylab="", xlim=c(117,176), main="2018", ylim=c(0,2.7), xaxt="n")
text(125,2.4, "n = 6791")
axis(side = 1, at = c(120,130,140,150,160,170), labels = c( '120','','','','','170'))
plot(as.numeric(as.character(sitemeans19$`sitemeans19$date`)),sitemeans19$`sitemeans19$cater`, pch=20, xlab="", ylab="", xlim=c(117,176), main="2019", ylim=c(0,2.7), xaxt="n")
text(125,2.4, "n = 8310")
axis(side = 1, at = c(120,130,140,150,160,170), labels = c( '120','','','','','170'))
plot(as.numeric(as.character(sitemeans20$`sitemeans20$date`)),sitemeans20$`sitemeans20$cater`, pch=20, xlab="Ordinal Date", ylab="", xlim=c(117,176), main="2020", ylim=c(0,2.7), xaxt="n")
text(125,2.4, "n = 3148")
axis(side = 1, at = c(120,130,140,150,160,170), labels = c( '120','','','','','170'))
plot(as.numeric(as.character(sitemeans$`sitemeans$date`)),sitemeans$`sitemeans$cater`, pch=20, xlab="Ordinal Date", ylab="", xlim=c(117,176), ylim=c(0,0.6), xaxt="n", yaxt="n")
text(126,2.4, "n = 31215")
axis(side = 1, at = c(120,130,140,150,160,170), labels = c( '120','130','140','150','160','170'))
axis(side = 2,
     ## Rotate the labels.
     las = 2)# saved as 10.7"x7.5"

par(mfrow=c(1,6),mar=c(4.1, 4.0, 2.5, 1), cex=0.9, cex.lab=1.25, cex.axis=1.2)
plot(as.numeric(as.character(sitemeans14$`sitemeans14$date`)),sitemeans14$`sitemeans14$cater`, pch=20, xlab="", ylab="Mean caterpillar abundance", xlim=c(117,176), main="2014", ylim=c(0,2.7), xaxt="n")
axis(side = 1, at = c(120,145,170), labels = c( '120','','170'))
text(145,2.7, "n=1106")
plot(as.numeric(as.character(sitemeans15$`sitemeans15$date`)),sitemeans15$`sitemeans15$cater`, pch=20, xlab="", ylab="", xlim=c(117,176), main="2015", ylim=c(0,2.7), xaxt="n")
axis(side = 1, at = c(120,145,170), labels = c( '120','','170'))
text(145,2.7, "n = 2632")
plot(as.numeric(as.character(sitemeans16$`sitemeans16$date`)),sitemeans16$`sitemeans16$cater`, pch=20, xlab="Ordinal Date", ylab="", xlim=c(117,176), main="2016", ylim=c(0,2.7), xaxt="n")
axis(side = 1, at = c(120,145,170), labels = c( '120','','170'))
text(145,2.7, "n = 2351")
plot(as.numeric(as.character(sitemeans17$`sitemeans17$date`)),sitemeans17$`sitemeans17$cater`, pch=20, xlab="", ylab="", xlim=c(117,176), main="2017", ylim=c(0,2.7), xaxt="n")
axis(side = 1, at = c(120,145,170), labels = c( '120','','170'))
text(145,2.7, "n = 6877")
plot(as.numeric(as.character(sitemeans18$`sitemeans18$date`)),sitemeans18$`sitemeans18$cater`, pch=20, xlab="", ylab="", xlim=c(117,176), main="2018", ylim=c(0,2.7), xaxt="n")
axis(side = 1, at = c(120,145,170), labels = c( '120','','170'))
text(145,2.7, "n = 6791")
plot(as.numeric(as.character(sitemeans19$`sitemeans19$date`)),sitemeans19$`sitemeans19$cater`, pch=20, xlab="", ylab="", xlim=c(117,176), main="2019", ylim=c(0,2.7), xaxt="n")
axis(side = 1, at = c(120,145,170), labels = c( '120','','170'))
text(145,2.7, "n = 8310") # saved as 5"x12"
