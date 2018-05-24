setwd('/Users/Nicola/Dropbox/EPFL_MAsem4/1_ATomistic_SIMulation/ATSIM/LAB4/results/1C')

gg <- read.table(file = "msd.dat",col.names = c("time","VAF","ddt"), comment.char = "#")
plot(gg$time,gg$VAF)
stdev <- c()
for(i in 1:length(gg$time)){
  stdev <- rbind(stdev,sd(gg$VAF[1:i]))
}
gg$sd <- stdev
plot(gg$time,gg$sd)
lines(c(0,5),c(1,1))
# useless -----------------------------------------------------------------

filenames = Sys.glob("*.dat")
# filename=filenames[1]
for(filename in filenames){
  gg <- read.table(file = filename,col.names = c("t","T","U","K"), comment.char = "#")
  f.spl = strsplit(filename,"_")[[1]]
  size=f.spl[4]
  step=strsplit(f.spl[6],".dat")[[1]]
  # cat("Size: ",size,"x",size,"x",size,"    Range:",range(gg$t),"\n")
  plot(gg$t,gg$T,
       main = paste0("Size: ",size,"x",size,"x",size,"   Step = ",step,"ps"),
       sub=paste0("Standard deviation:",sd(gg$T)/mean(gg$T)*100," %"),
       ylab = "Total energy",xlab="time step",
       type="l",col="red")
}
