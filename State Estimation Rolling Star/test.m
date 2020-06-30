clc
clear all
close all

syms t1 t2 a1 a2

[T0_n,Tnm1_n] = fwdkinRISS([a1,a2], [0,0], [0,0], [t1,t2]);

P01 = T0_n{1}(1:3,4)
x20 = T0_n{2}(1:3,1)
P02g = P01+(a2/2)*x20
p1g2