#!/usr/bin/env julia

# about the problem
mdMax = 10000
m = 3.0
k = 2.0
dt = 0.001

# initial conditions

x = -1.
v = 0.
a = 0.
t = 0.


for i = 1:mdMax
    f = -k*x
    a = f/m
    v += a*dt
    x += v*dt
    t += i*dt

    println("$t $x")
end
