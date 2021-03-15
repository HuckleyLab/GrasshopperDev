#Q10: growth~1, development ~2

#long day: 14/10
#short day: 12/12

#LV: 24 +- 2 
#HV: 24+- 4

#test out Q10 equation
#Q10=(R2/R1)^(10/(T2-T1))
R2= function(T2, R1=1, Q10=1, T1=1) R1 *Q10^((T2-T1)/10)

#plot out rates with different Q10s
r2s=sapply(1:50, FUN="R2", Q10=1.5,T1=24)
plot(1:50, r2s, type="l")
r2s=sapply(1:50, FUN="R2", Q10=2,T1=24)
points(1:50, r2s, type="l", col="red")

#---
#functions for days development and growth
#assumes R1=1 and arbitrarily requires 1000 units of development
#12 hours of each temperature (T2l and T2h)
days.dev=function(Q10, T2l=22, T2h=26, T1=22)  1/(12/1000*(Q10^((T2l-T1)/10)+ Q10^((T2h-T1)/10)))
growth=  function(Q10, days, T2l=22, T2h=26, T1=22)  12*(Q10^((T2l-T1)/10)+ Q10^((T2h-T1)/10))*days

#plot out potential development times across Q10s for lv and hv scenario
#assume T1=24 to produce observed dynamics
days.hv=sapply(seq(0.6,3,0.1), FUN="days.dev", T2l=20, T2h=28, T1=24) 
plot(seq(0.6,3,0.1), days.hv, type="l", col="red")
days.lv=sapply(seq(0.6,3,0.1), FUN="days.dev", T2l=22,T2h=26, T1=24)
points(seq(0.6,3,0.1), days.lv, type="l", col="blue")

#plot out potential growth across Q10s for lv and hv scenario
growth.hv= sapply(1:5, FUN="growth", days= days.dev(Q10=2, T2l=20, T2h=28, T1=24),T2l=20, T2h=28, T1=24) 
plot(1:5, growth.hv, type="l", col="red")
growth.lv= sapply(1:5, FUN="growth", days= days.dev(Q10=2, T2l=22, T2h=26, T1=24), T2l=22, T2h=26, T1=24) 
points(1:5, growth.lv, type="l", col="blue")

#Plot out example Q10 scenarios for A1 and C1 sites
#development
days.hv.a1= days.dev(Q10=1.5, T2l=20, T2h=28, T1=24)
days.hv.c1= days.dev(Q10=2, T2l=20, T2h=28, T1=24)
days.lv.a1= days.dev(Q10=1.5, T2l=22, T2h=26, T1=24)
days.lv.c1= days.dev(Q10=2, T2l=22, T2h=26, T1=24)

plot(c(2195,3048), c(days.hv.a1,days.hv.c1), type="l", col="red", ylim=c(30,45))
points(c(2195,3048), c(days.lv.a1,days.lv.c1), type="l", col="blue")

#growth
growth.hv.a1= growth(Q10=3, days=days.hv.a1, T2l=20, T2h=28, T1=24)
growth.hv.c1= growth(Q10=2.5, days=days.hv.c1, T2l=20, T2h=28, T1=24)
growth.lv.a1= growth(Q10=3, days=days.lv.a1, T2l=22, T2h=26, T1=24)
growth.lv.c1= growth(Q10=2.5, days=days.lv.c1, T2l=22, T2h=26, T1=24)

plot(c(2195,3048), c(growth.hv.a1,growth.hv.c1), type="l", col="red", ylim=c(1000,1100))
points(c(2195,3048), c(growth.lv.a1,growth.lv.c1), type="l", col="blue")








