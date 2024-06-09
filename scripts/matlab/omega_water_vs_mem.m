clc
clear all
close all

g = 9.81;

h = 10;
m_by_rho = 10*0.9;
Ten_by_rho = 5*g*10;

L = 2/4*20;
k = 2*pi/L;

ww = sqrt(g*k*tanh(k*h));

syms wm
wm=solve(wm*wm - g*k*tanh(k*h)*(1 - m_by_rho/g*wm*wm + Ten_by_rho/g*k*k));
wm = abs(double(wm));

fprintf('Lw = %f \n', L)
fprintf('ww = %f \n', ww)
fprintf('\n')
fprintf('Lm = %f \n', L)
fprintf('wm = %f \n', wm)