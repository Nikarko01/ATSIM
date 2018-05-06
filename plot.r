mgo <- read.table("~/Dropbox/EPFL_MAsem4/1_ATomistic_SIMulation/ATSIM/pb1mgloop3/Mg_results.dat")


library(plotly)

plot_ly(data = mgo,x=~V1,y=~V2,z=~V3)

minmgo <- subset(mgo,V3==min(mgo$V3))
minV1 <- subset(mgo,V1==1.61)

plot(minV1$V2,minV1$V3)
c <- polyfit(minV1$V2,minV1$V3, 2)
x=minV1$V2
x=seq(0,10,length.out = 16)
y <- c[1]*x^2+c[2]*x+c[3]
plot(x,y,col="green")
library(pracma)

-c[2]/2/c[1]

# V1   V2       V3
# 30 1.61 6.05 -250.648  

