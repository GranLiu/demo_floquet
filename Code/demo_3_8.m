% -------------------------------------------------------------------------
% Realize Fig.3.8 in Phased Array Antennas.
% 
% Yongxi Liu, Xi'an Jiaotong University, 2023-03.
% -------------------------------------------------------------------------
clc;
clear;
close all;

%% simulation params
f = 2e9;
lbd = 3e8/f;
k0 = 2*pi/lbd;
a = (0.5:0.01:1.5)*lbd;
b = 0.2*lbd;
h = 0.05*lbd;
ep_r = 2.0;

krou_init = k0*(1+0.5*(k0*h*(1-1/ep_r))^2);
syms x
eqn = sqrt(k0^2-x^2)+1j*sqrt(ep_r*k0^2-x^2)*tan( sqrt(ep_r*k0^2-x^2)*h )/ep_r == 0;
S = vpasolve(eqn,x,krou_init);
krou = double(S);

theta_ep1 = zeros(size(a));
for i=1:length(a)
    if(2*pi/a(i) > krou)
        theta_ep1(i) = asin( (2*pi/a(i) - krou)/k0 )/pi*180;
    else
        theta_ep1(i) = asin( (krou - 2*pi/a(i))/k0 )/pi*180;
    end
end

ep_r = 3.0;
krou_init = k0*(1+0.5*(k0*h*(1-1/ep_r))^2);
syms x
eqn = sqrt(k0^2-x^2)+1j*sqrt(ep_r*k0^2-x^2)*tan( sqrt(ep_r*k0^2-x^2)*h )/ep_r == 0;
S = vpasolve(eqn,x,krou_init);
krou = double(S);

theta_ep2 = zeros(size(a));
for i=1:length(a)
    if(2*pi/a(i) > krou)
        theta_ep2(i) = asin( (2*pi/a(i) - krou)/k0 )/pi*180;
    else
        theta_ep2(i) = asin( (krou - 2*pi/a(i))/k0 )/pi*180;
    end
end
% 

figure();hold on;
plot(a/lbd,theta_ep1,'k');
plot(a/lbd,theta_ep2,'--b');
legend('$\epsilon_\mathrm{r}=2.0$','$\epsilon_\mathrm{r}=3.0$','interpreter','latex','fontsize',12);
xlabel('$a/\lambda$','interpreter','latex','fontsize',12);
ylabel('Blind angle [deg]','interpreter','latex','fontsize',12);
ylim([0,90])
xlim([0.5,1.5])
grid on;
set(gca,'LooseInset',get(gca,'TightInset'));

figure(); hold on;
ang = 0:0.01:2*pi;
plot(krou*cos(ang),krou*sin(ang));
plot(k0*cos(ang),k0*sin(ang));
grid on;
axis equal;

