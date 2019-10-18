import numpy
from numpy import linalg as LA
from sympy import Matrix
import math

k = numpy.array([[0.944444, 0.5],[0.0555555, 0.5]]) 
print(k)
w, v = LA.eig(k)
print(w)
print(v)

k2 = numpy.array([[0.9, 0.9],[0.1, 0.1]]) 
print(k2)
w, v = LA.eig(k2)
print(w)
print(v)
test = LA.matrix_power(k2, 10); 
print(test)

x0 = numpy.array([0.5,0.5]);
x = x0;
for i in range(10):
    print(x)
    x = numpy.dot(k,x);

X1 = numpy.array(
        [[0.666667,  0.333333,         0,  0.333333], 
         [0.333333,  0.333333, 0.0833333,         0], 
         [0,  0.333333,      0.75,  0.333333], 
         [0,         0,  0.166667,  0.333333]]) 
print(X1)
X1 = Matrix(X1)

P, J = X1.jordan_form()
print("P")
print(P)

print("J")
print(J)

#X1 = numpy.array([0.9, 0.05, 0.05]);
#K1 = numpy.array(
#    [[0.9, 0.9, 0.90], 
#    [0.05, 0.05, 0.05], 
#    [0.05, 0.05, 0.05]])
#
#print(K1)
#x = X1
#for i in range(10):
#    print(x)
#    x = numpy.dot(K1,x);
