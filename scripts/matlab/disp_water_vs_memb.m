clc
clear all
close all


h = 10;
w = 0.5:0.05:3.5;
T = 2*pi./w;

kw = zeros(size(w));
km = zeros(size(w));


for i = 1:size(w,2)
    kw(i) = disp_water(w(i), h);
    km(i) = disp_memb(w(i), h);
end

Lw = 2*pi./kw;
Lm = 2*pi./km;

save("res.mat")

figure(1)
hold on
plot(w, kw, 'k', 'LineWidth', 3)
plot(w, km, 'r', 'LineWidth', 3)
grid on
xlabel('w (rad/s)')
ylabel('k (rad/m)')
legend('Water','Membrane')
set(gca,'GridAlpha',1,'GridLineStyle','--')
xline(0.7,'HandleVisibility','off')


figure(2)
hold on
plot(T, Lw, 'k', 'LineWidth', 3)
plot(T, Lm, 'r', 'LineWidth', 3)
grid on
xlabel('T (s)')
ylabel('L (m)')
legend('Water','Membrane')
set(gca,'GridAlpha',1,'GridLineStyle','--')
xline(2*pi/0.7,'HandleVisibility','off')


figure(3)
hold on
plot(w, w./kw, 'k', 'LineWidth', 3)
plot(w, w./km, 'r', 'LineWidth', 3)
grid on
xlabel('w (rad/s)')
ylabel('c (m/s)')
legend('Water','Membrane')
set(gca,'GridAlpha',1,'GridLineStyle','--')
xline(0.7,'HandleVisibility','off')


%% Functions
function k = disp_water(w, h)
    g=9.81;
    T=2*pi/w;
    
    syms L
    L=solve(L-(g/2/pi*T*T*tanh(2*pi/L*h)));
    L=abs(double(L));
    
    k = 2*pi/L;
end


function k = disp_memb(w, h)
    g=9.81;
    m_by_rho = 0.9;
    Ten_by_rho = g*10;
    
    syms k
    k=solve(w*w - g*k*tanh(k*h)*(1 - m_by_rho/g*w*w + Ten_by_rho/g*k*k));
    k=abs(double(k));        
    k = k(1);
    
end
