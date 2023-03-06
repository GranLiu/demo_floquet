% -------------------------------------------------------------------------
% Realize Fig.3.4 & Fig.3.5 in Phased Array Antennas.
% 
% Yongxi Liu, Xi'an Jiaotong University, 2023-03.
% -------------------------------------------------------------------------
clc;
% clear;
close all;

%% simulation params
f = 2e9;
lbd = 3e8/f;
k0 = 2*pi/lbd;
a = (0.5:0.01:1.5)*lbd;
b = 0.2*lbd;
h = 0.05*lbd;
ep_r = 2.0;

% krou = k0*(1+0.5*(k0*h*(1-1/ep_r))^2);


syms x
eqn = sqrt(k0^2-x^2)+1j*sqrt(ep_r*k0^2-x^2)*tan( sqrt(ep_r*k0^2-x^2)*h )/ep_r == 0;
S = vpasolve(eqn,x,k0);
krou = double(S);

theta = asin( (2*pi./a - krou)/k0 )/pi*180;
theta2 = asin( (2*2*pi)./a - krou)/k0 )/pi*180;

figure();hold on;
plot(a/lbd,theta);
plot(a/lbd,theta2);
grid on;

ep_r = 3.0;
syms x
eqn = sqrt(k0^2-x^2)+1j*sqrt(ep_r*k0^2-x^2)*tan( sqrt(ep_r*k0^2-x^2)*h )/ep_r == 0;
S = vpasolve(eqn,x,k0);
krou = double(S);


% krou = k0*(1+0.5*(k0*h*(1-1/ep_r))^2);
theta = asin( (2*pi./a - krou)/k0 )/pi*180;
% plot(a/lbd,theta);

ylim([0,90])

figure(); hold on;
ang = 0:0.01:2*pi;
plot(krou*cos(ang),krou*sin(ang));
plot(k0*cos(ang),k0*sin(ang));
grid on;
axis equal;
%% E-plane scan

% sqrt(k0^2-krou^2)+1j*sqrt(ep_r*k0^2-krou^2)*tan( sqrt(ep_r*k0^2-krou^2)*h )/ep_r
% sqrt(k0^2-zz^2)+1j*sqrt(ep_r*k0^2-zz^2)*tan( sqrt(ep_r*k0^2-zz^2)*h )/ep_r
