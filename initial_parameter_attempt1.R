library(deSolve)
library(phaseR)
library(pracma)
RHS_op <- function(t,y,parameters){
  G <- y[1]  # the first state variable is E, the effector cells
  U <- y[2]   # and the second is T, the tumor cells
  S <- y[3]
  F <- y[4]
  T <- y[5]
  R <- y[6]
  a <- parameters[1]  # this function has 8 parameters
  b <- parameters[2]
  c <- parameters[3]
  d <- parameters[4]
  e <- parameters[5]
  f <- parameters[6]
  g <- parameters[7]
  d1 <- parameters[8]
  e1 <- parameters[9]
  h <- parameters[10]
  dy <- numeric(6)    # the output will be a numeric vector of length 6, because we have 6 differential equations
  dy[1] <- -a*G*U
  dy[2] <- a*G*U + b*R*U - c*U - d*U - e*U + f*T - g*U*S
  dy[3] <- g*U*S - d1*S - e1*S
  dy[4] <- c*U
  dy[5] <- d*U - h*T - f*T + d1*S
  dy[6] <- h*T + e*U - b*R*U + e1*S
  list(dy)  # this tells R what to output: the vector of derivatives
}

#now I will set values of parameters
a = 0.09
b = (0.01+.2)/2 #i took the avg of the two relapse parameters from zoe's notes
c =  0.65 #https://link.springer.com/article/10.1007/s11538-022-01002-w/tables/2 (took most recent value from 2019)
d = 0.4 #from article saying around of 45% of users experience non fatal overdose linked in overleaf
e = (0.01+.2)/2 #average of two recovery rates from zoe's notes
f = 0.85 #average of recovery rate from zoe's notes
g = 0.15 #from mira's estiation
d1 = 0.001 #one non-fatal per 1000 injections, makes sense if we assume one injection per perso, looking for information about injections/person
e1 = 0.2 #ill say in an SIS, double the % of people would choose recovery
h = 0.25 #from statistic Sam said

#Now I need to solve the DEs - trying to use dr. rad's tumor model code

systemsoln = function(parameters,x0,tn){
  parameters = c(a,b,c,d,e,f,g,d1,e1,h)
  tn = seq(0, 100, by = 0.01)
  x0 = c(G = 2000,U = 1000,S = 1,F = 0,T = 1,R = 1)
  Y = ode(func = RHS_op, times = tn, y = x0, parms = pars);
  return(Y)
}


# going to try a slightly different approach using code from online

parameters <- c(a = 0.09,b=0.1,c=0.0283,d=0.4,e=0.1,f=0.75,g=0.15,d1=0.2,e1=0.2,h=0.25)
state <- c(G = 2000,U = 1000,S = 1,F = 0,T = 1,R = 1)

RHS_op <- function(t,state,parameters){
  with(as.list(c(state,parameters)),{
    dG <- -a*G*U
    dU <- a*G*U + b*R*U - c*U - d*U - e*U + f*T - g*U*S
    dS <- g*U*S - d1*S - e1*S
    dF <- c*U
    dT <- d*U - h*T - f*T + d1*S
    dR <- h*T + e*U - b*R*U + e1*S
    list(c(dG,dU,dS,dF,dT,dR))
  })
}

times <- seq(0,100,by = 0.01)

Y <- ode(y = state, times = times, func = RHS_op, parms = parameters)
head(Y) 

#this seems to have given something over time! so now we can try to plot the results and see what happens

par(oma = c(0, 0, 3, 0))
plot(Y, xlab = "time", ylab = "-")
mtext(outer = TRUE, side = 3, "Opioid Model", cex = 1.5)

