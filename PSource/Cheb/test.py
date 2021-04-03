"""
Simple routines to test cheb module.
"""

import cheb
import numpy as np

def norm(arr):
   return sum(abs(arr))/len(arr)

def der_diff(n,lower,upper):
   def f(x):
      return 1.0 / (1.0 + pow(x,2))
   def df(x):
      return - 2.0 * x / pow(1.0 + pow(x,2),2)

   cheb.init(n,lower,upper)
   v   = np.array([f(cheb.pt(i))  for i in range(cheb.n())])
   dv  = np.array([0.0            for i in range(cheb.n())])
   dv2 = np.array([df(cheb.pt(i)) for i in range(cheb.n())])
   cheb.der(v,dv)
   print(norm(dv-dv2))
   cheb.cleanup()

def filter_change(n,lower,upper):
   import random as rand

   cheb.init(n,lower,upper)
   v = np.array([0.1*(-0.5+rand.random()) for i in range(cheb.n())])
   print(norm(v))
   cheb.filter(v)
   print(norm(v))
   cheb.cleanup()


print("Test 1: each should be decreasing")
der_diff(32 ,-1.1,52.7)
der_diff(64 ,-1.1,52.7)
der_diff(128,-1.1,52.7)
der_diff(256,-1.1,52.7)
print("Test 2: each should be decreasing")
der_diff(32 ,52.7,0.0)
der_diff(64 ,52.7,0.0)
der_diff(128,52.7,0.0)
der_diff(256,52.7,0.0)
print("Test 3: each should be decreasing")
filter_change(32 ,52.7,0.0)
print("---")
filter_change(64 ,52.7,0.0)
print("---")
filter_change(128,52.7,0.0)
print("---")
filter_change(256,52.7,0.0)
