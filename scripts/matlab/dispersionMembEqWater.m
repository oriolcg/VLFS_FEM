clc
clear all 
close all

g=9.81;

m_by_rho = 10*0.9;
Ten_by_rho = 5*g*10;

%% k from w
h=10;
syms k
k=solve(k-m_by_rho/Ten_by_rho*g*tanh(k*h));
k=abs(double(k));

fprintf('m/Ten = %f\n',m_by_rho/Ten_by_rho);
fprintf('m/Ten*g = %f\n',m_by_rho/Ten_by_rho*g);
fprintf('L > %f\n',Ten_by_rho/m_by_rho/g*2*pi);
fprintf('k < %f\n',m_by_rho/Ten_by_rho*g);
fprintf('k = %f\n',k);
fprintf('L = %f\n',2*pi/k);


k=0.17
k-m_by_rho/Ten_by_rho*g*tanh(k*h)

