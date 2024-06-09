
clc
clear all 
close all

g=9.81;

T=2*pi/2.51;
h=10;%0.4572; %0.1524
amp=0.0875/2;
syms L
L=solve(L-(g/2/pi*T*T*tanh(2*pi/L*h)));
L=abs(double(L));

hByL=h/L;
fprintf('f = %f\n',1/T);
fprintf('T = %f\n',T);
fprintf('Amp = %f\n',amp);
fprintf('L = %f\n',L);
fprintf('h = %f\n',h);
fprintf('h/L = %f\n',hByL);
fprintf('h/L0 = %f\n',h/(g/2/pi*T*T));
fprintf('kh = %f\n',2*pi*hByL);
fprintf('w = %f\n',2*pi/T);
fprintf('UMax = %f\n',amp*g*T/L);
fprintf('UMax x T = %f\n',amp*g*T/L*T);
fprintf('h/g/T/T = %f\n',h/g/T/T);
fprintf('H/g/T/T = %f\n',2*amp/g/T/T);  

k=2*pi/L;
kh=k*h;
tfPis=2*(cosh(2*kh)-1)/(sinh(2*kh) + 2*kh);
tfPis=1/tfPis;
fprintf('S/H = %f\n',tfPis);
fprintf('s = %f\n',amp*tfPis);

%% Scaled
sc=8.0392;
scsq=sqrt(sc);
T=T*scsq;
amp=amp*sc;
h=h*sc;
LEst=L*sc;

syms L
L=solve(L-(g/2/pi*T*T*tanh(2*pi/L*h)));
L=abs(double(L));

hByL=h/L;
fprintf('f = %f\n',1/T);
fprintf('T = %f\n',T);
fprintf('Amp = %f\n',amp);
fprintf('L Sca = %f\n',LEst);
fprintf('L Cal = %f\n',L);
fprintf('h = %f\n',h);
fprintf('h/L = %f\n',hByL);
fprintf('UMax = %f\n',amp*g*T/L);
fprintf('UMax x T = %f\n',amp*g*T/L*T);