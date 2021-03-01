#### playing with normal distributions ####

# y= (1/(σ*sqrt(2*pi)))*exp((-(x-μ)^2)/(2*σ^2))   equation for normal distribution

x=seq(-50,50,0.5)
σ=4
μ=2
h=0.1
y <- (1/(σ*sqrt(2*pi)))*exp(((-(x-μ)^2)/(2*σ^2))) 
plot(x,y, type="l")

x=seq(-50,50,0.5)
σ=3
μ=2
h=0.1
y <- exp(((-(x-μ)^2)/(2*σ^2))+log(h))
y2 <- exp(((-(x-μ)^2)/(2*σ^2))+log(h))
plot(x,y, type="l", ylim=c(0,0.5))
points(x,y2, type="l")

# (1/(σ*sqrt(2*pi))) I think this part makes the area underneath the curve = 1

## exp(((-(x-μ)^2)/(2*σ^2))) 
#change in mu =direct shift if peak position, width and height unaffected
#change in sigma =direct shift in width, position and height unaffected, sign irrelevant
#when h is introduced an increase in h will increase width, needs to be positive and more than 0
#needs to be log(h) for direct shift in height (width still increases)  

# we want mu, sigma and h to interact with temperature

y <- caterpillar abundance
x <- date
z <- temperature
μ <- peak date (a1+b1*z+e1) 
σ <- standard deviation(a2+b2*z+e2)
h <- height (has to be >0)(a3+b3*z+e3)

y = exp(( (-(x-(a+b1*z))^2) / (2*σ^2) ) +log(h))



