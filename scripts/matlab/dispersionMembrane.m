clc
clear all 
close all

g=9.81;

m_by_rho = 0.9;
Ten_by_rho = g*10;

%% k from w
w = 2.51;
T = 2*pi/w;
h=10;
amp=0.01;
syms k
k=solve(w*w - g*k*tanh(k*h)*(1 - m_by_rho/g*w*w + Ten_by_rho/g*k*k));
k=abs(double(k));

fprintf('w = %f\n',w);
fprintf('f = %f\n',1/T);
fprintf('T = %f\n',T);
fprintf('k = %f\n',k);
fprintf('L = %f\n',2*pi/k);
fprintf('h = %f\n',h);
fprintf('kh = %f\n',k*h);
fprintf('Dispersion Factor = %f \n',(1 - m_by_rho/g*w*w + Ten_by_rho/g*k*k))


%% w from k

% k = 2*pi/20;
% h=10;
% amp=0.01;
% syms w
% w=solve(w*w - g*k*tanh(k*h)*(1 - m_by_rho/g*w*w + Ten_by_rho/g*k*k));
% w=abs(double(w));
% T = 2*pi./w;
% 
% fprintf('w = %f\n',w);
% fprintf('f = %f\n',1./T);
% fprintf('T = %f\n',T);
% fprintf('k = %f\n',k);
% fprintf('L = %f\n',2*pi/k);
% fprintf('h = %f\n',h);
% fprintf('kh = %f\n',k*h);
% fprintf('Dispersion Factor = %f \n',(1 - m_by_rho/g*w.*w + Ten_by_rho/g*k*k))
