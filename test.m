clear all
clc
%mex -largeArrayDims -O -c smul.c
mex -largeArrayDims smul.c
A=sprand(100,100,0.001);
B=sprand(100,100,0.001);
C=sadd(A,B);
%(C - (A+B))
