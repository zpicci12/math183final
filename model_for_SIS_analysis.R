library(deSolve)
library(phaseR)
library(pracma)

#these parameters maintain a slight growth in S, 
#2% of Susceptible are Abusers
#The recovery population is 1/2 the Abuser population

a = 0.00000015 # from susceptible to abuse/misuse
b = 0.0000000025 # interaction term between susceptible and abuse
c = 0.38 # from abuse to recovery
d = 0.008 # death of abusers/misusers
e = 0.75 #relapse
g = 0.4 #from SIS abuse to recovery (assuming people at SIS recover slightly more than general abusers)
l = 3  #adjusting this parameter would adjust willingness to use an SIS and determine how quickly U reaches carrying capacity
k = 6000 #adjusting this parameter changes the SIS carrying capacity

##Assumptions for LA County Population: 
#S0: 9,830,000 (2021 LA County Population) #find adult population? 
#A0: 196,600 (estimate that 2% of adults abuse opioids)
#parameter fitting: find parameters such that given S0 and A0, there are D1 deaths 
#There are, on average 1,570 fatal overdoses per year (without SIS)
#Assume at time t=0, 100 people use SIS

state <- c(S = 9830000, A = 196600, R = 98300, D = 0, U = 100)
parameters <- c(a,b,c,d,e)

N = 9830000

pre_op <- function(t,state,parameters){
  with(as.list(c(state,parameters)),{
    dS <- -(b*A)*(S) - b*S*(A+U/2) + (0.00055*N) #added a natural replacement rate for better long-term outcomes
    dA <- b*S*(A+U/2) - (c+d)*A + e*R + (b*A)*S - 0.0118*A - l*U*(1-(U/k)) 
    dR <- c*A - e*R - 0.00372*R + g*U 
    dD <- d*A
    dU <- l*U*(1-(U/k)) - g*U 
    list(c(dS,dA,dR,dD, dU))
  })
}

times <- seq(0,30, 0.001)

Y <- ode(y = state, times = times, func = pre_op, parms = parameters)
head(Y) 
tail(Y)

par(oma = c(0, 0, 3, 0))
plot(Y, xlab = "time", ylab = "-", xlim=c(0,30))
mtext(outer = TRUE, side = 3, "Opioid Model", cex = 1.5)
