rm(list=ls())
#
library(reticulate)     # for reshaping
#
# io data info
# [, 1: 4]  k,pc(k,nt),tc(k,nt),wc(k,nt), &
# [, 5: 6]  qv(i,k),cldn(k,nt),  &  
# [, 7:11]  qc(k,nt),qi(k,nt),qs(k,nt),qr(k,nt),qg(k,nt), &
# [,12:16]  nc(k,nt),ni(k,nt),ns(k,nt),nr(k,nt),ng(k,nt)  # <-- input ends
# [,17:18]  effcs(i,k,nt),effis(i,k,nt)
#
kMax=42     # vertical levels
rlev=kMax   # required levels
deltat=450.0  # in seconds
#
# 24*3600/deltat = 192 timesteps per day
# 2304/192 = 12 days
# 12*192*kMax = 96768 rows [42 vertical levels]
# 12*192*42-97608 = -840 rows
#
inp<-read.table('invars_T126_L042_lat10.70S_lon291.dat',header = F,skip=840)
out1<-read.table('outvars_iact1_T126_L042_lat10.70S_lon291.dat',header = F,skip=840)
out2<-read.table('outvars_iact2_T126_L042_lat10.70S_lon291.dat',header = F,skip=840)
#
# number mixing ratio & reshaping to [n3hrs x ndeltatper3hrs x vlevs]
nc0=array_reshape(inp[,12],c(12*8,24,42))
nc1=array_reshape(out1[,12],c(12*8,24,42))
nc2=array_reshape(out2[,12],c(12*8,24,42))
cld=array_reshape(inp[,6],c(12*8,24,42))
#
# calculate 3 hourly mean
nc0_3hrly=apply(nc0,c(1,3), mean.default)
nc1_3hrly=apply(nc1,c(1,3), mean.default)
nc2_3hrly=apply(nc2,c(1,3), mean.default)
cld_3hrly=apply(cld,c(1,3), mean.default)
#
plot(seq(1,96),nc0_3hrly[,18],type="b",xlab='num. 03 hrly steps for 12 days',ylab='number mixing ratio (NC)',ylim=c(0,max(nc1_3hrly)),xlim=c(0,97))
lines(seq(1,96),nc2_3hrly[,18],type="b",col="green")
#lines(seq(1,96),nc1_3hrly[,18],type="b",col="red")
#
#



