#!/usr/bin/env python

mdMax = 10000
# about the problem
m = 3.0
k = 2.0
dt = 0.001

# initial conditions

i = 0
x = -1. 
v = 0.
a = 0.
t = 0.


while i < mdMax:
    f = -k*x
    a = f/m
    v += a*dt
    x += v*dt
    t += i*dt
    i += 1

#    print t, x
print 'done'
