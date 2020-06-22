%%RISS 2020 STATICS 4-BAR WHOLE%%

clc
clear all
close all 

syms f1 f2 a1 a2 a3 a1p t1p m1 m2 m3 m1p g x2 y2

eq1 = 2*a1p*cos(t1p)*x2+2*a1p*sin(t1p)*y2 == a1p^2-a3^2+x2^2+y2^2
eq2 = a2^2 == a1^2-2*a1*x2+y2^2

[solx,soly] = solve(eq1, eq2)