#########################################
## Find the optimum t in c(0, T) for   ## 
## 	h(t;q,k) = (kq)^t + q^{T-t}-k^tq^T ## 
#########################################

#install.packages("ggplot2")
library(ggplot2)

ht <- function(t,q,k,Tx) (k*q)^t*log(k*q) - q^(Tx-t)*log(q)-k^t*q^Tx*log(k)
htprime <- function(t,q,k,Tx) {(k*q)^t*log(k*q)^2 + (q^(Tx-t))*log(q)^2 -(k^t)*(q^Tx)*log(k)^2 }


NR.cap <- function(f, fprime, a, b, tol = 1e-8, n = 10000,...) {
  
  x0 <- a # Set start value to supplied lower bound
  d <- n # Initialize for iteration results
  
  # Check the upper and lower bounds to see if approximations result in 0
  fa <- f(a,...)
  if (fa == 0.0) {
    return(a)
  }
  
  fb <- f(b,...)
  if (fb == 0.0) {
    return(b)
  }
  
  for (i in 1:n) {
    x1 <- x0 - (f(x0,...) / fprime(x0,...)) # Calculate next value x1
    d[i] <- x1 # Store x1
    # Once the difference between x0 and x1 becomes sufficiently small, output the results.
    if (abs(x1 - x0) < tol) {
      t.opt <- tail(d, n=1)
      #result <- list('t.optimum' = t.opt, 'iterations' = d)
      return(t.opt)
    }
    
    x0 <- x1
  }
  # If Newton-Raphson has not yet reached convergence set x1 as x0 and continue 
  print('Too many iterations in method')
}


numq=seq(0.1, 0.9, by = 0.1) 
numk=seq(0.1, 2.0, by = 0.1)
TT=10
xx <- matrix(nrow = length(numk),ncol=length(numq))
for(i in 1:length(numk)){
  for(j in 1:length(numq)){
    kk<-numk[i]
    qq <- numq[j]
    if(kk*qq>0 & kk*qq <1){
      topt <- NR.cap(ht,htprime,a=TT/2,b=TT,q=qq,k=kk,Tx=TT)
      xx[i,j] <- topt}
  }}

##################
##Plot t optimum##
##################

TT <- 10
yt <- data.frame(k=numk,xx)
names(yt)<-c("k",paste0("q",numq))
print(yt)

colors <- c("q=0.1"="black", "q=0.2"="red", "q=0.3"="green", 
            "q=0.4"="blue", "q=0.5"="purple", "q=0.6"="brown",
            "q=0.7"="blue4", "q=0.8"="orange", "q=0.9"="deeppink")

opt.yt_plot <- ggplot(yt, aes(k)) +  
  geom_line(aes(y=q0.1, color = "q=0.1"),size = 0.5) +
  geom_line(aes(y=q0.2, color = "q=0.2"),size = 0.5) +
  geom_line(aes(y=q0.3, color = "q=0.3"),size = 0.5) +
  geom_line(aes(y=q0.4, color = "q=0.4"),size = 0.5) +
  geom_line(aes(y=q0.5, color = "q=0.5"),size = 0.5) +
  geom_line(aes(y=q0.6, color = "q=0.6"),size = 0.5) +
  geom_line(aes(y=q0.7, color = "q=0.7"),size = 0.5)  +
  geom_line(aes(y=q0.8, color = "q=0.8"),size = 0.5) +
  geom_line(aes(y=q0.9, color = "q=0.9"),size = 0.5) +
  labs(x="k", y="t optimum", color ="Legend") +
  scale_color_manual(values = colors) +
  theme(legend.position="bottom")

opt.yt_plot


##################
##Plot ratio t/T##
##################

TT <- 10
rt <- data.frame(k=numk,xx/TT)
names(rt)<-c("k",paste0("q",numq))
print(rt)

colors <- c("q=0.1"="black", "q=0.2"="red", "q=0.3"="green", 
            "q=0.4"="blue", "q=0.5"="purple", "q=0.6"="brown",
            "q=0.7"="blue4", "q=0.8"="orange", "q=0.9"="deeppink")

opt.rt_plot <- ggplot(rt, aes(k)) +  
  geom_line(aes(y=q0.1, color = "q=0.1"),size = 0.5) +
  geom_line(aes(y=q0.2, color = "q=0.2"),size = 0.5) +
  geom_line(aes(y=q0.3, color = "q=0.3"),size = 0.5) +
  geom_line(aes(y=q0.4, color = "q=0.4"),size = 0.5) +
  geom_line(aes(y=q0.5, color = "q=0.5"),size = 0.5) +
  geom_line(aes(y=q0.6, color = "q=0.6"),size = 0.5) +
  geom_line(aes(y=q0.7, color = "q=0.7"),size = 0.5)  +
  geom_line(aes(y=q0.8, color = "q=0.8"),size = 0.5) +
  geom_line(aes(y=q0.9, color = "q=0.9"),size = 0.5) +
  labs(x="k", y="t / T ratio", color ="Legend") +
  scale_color_manual(values = colors) +
  theme(legend.position="bottom")

opt.rt_plot
