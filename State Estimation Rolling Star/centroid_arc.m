%%RISS 2020 FWD KIN 4 LEGGED ROLLING STAR%%

clc
clear all
close all 

%knowns
L = 1;
k = 1:.1:(2*pi)/L;


alpha = (L*k)/2; 
z = (sin(alpha))./(k.*alpha) - cos(alpha)./k;

plot(k,z)
xlabel('k')
ylabel('dist to chord')

[M,I] = max(z)
k_max = k(I)
alpha_max = alpha(I)
