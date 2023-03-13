% -------------------------------------------------------------------------
% Realize Fig.3.9 in Phased Array Antennas.
% 
% Yongxi Liu, Xi'an Jiaotong University, 2023-03.
% -------------------------------------------------------------------------
clc;
clear;
close all;

addpath("./func");

%% simulation params
k0 = 2*pi;
a = 0.6;
b = 0.6;
l = 0.4;
t = 0.1;
h = 0.05;
ep_r = 2.0;
Z_0 = 377;
d_theta = 0.05;
theta = (d_theta:d_theta:90)/180*pi;

% source impedance is conjugate matched to bore-sight scan
[R_g,X_g] = z_fl_dp(0.01,0,ep_r,a,b,h,l,t);
X_g = -X_g;
I_inc = sqrt(4*pi/R_g);

% Embedded element gain for E-plane scan
G_e = zeros(size(theta));
for idx = 1:length(theta)
    theta_0 = theta(idx);
    phi_0 = 0;
    cos_t = cos(theta_0);
    sin_t = sin(theta_0);
    cos_p = cos(phi_0);
    sin_p = sin(phi_0);
    [R_fl,X_fl] = z_fl_dp(theta_0,phi_0,ep_r,a,b,h,l,t);
    nrt_p = sin_p*sinc( k0*(l/2)*sin_t*cos_p/pi ) *...
        sinc( k0*(t/2)*sin_t*sin_p/pi );
    dnrt_p = 1 - 1j*sqrt(ep_r-sin_t^2)*cot(k0*h*sqrt(ep_r-sin_t^2))/cos_t;
    E_p = 1j*I_inc*Z_0*l*(2*R_g/(R_fl+R_g+1j*(X_fl+X_g))) * nrt_p/dnrt_p;
    nrt_t = cos_t*cos_p*sinc( k0*(l/2)*sin_t*cos_p/pi ) *...
        sinc( k0*(t/2)*sin_t*sin_p/pi );
    dnrt_t = 1 - 1j*ep_r*cos_t/sqrt(ep_r-sin_t^2)*cot(k0*h*sqrt(ep_r-sin_t^2));
    E_t = -1j*I_inc*Z_0*l*(2*R_g/(R_fl+R_g+1j*(X_fl+X_g))) * nrt_t/dnrt_t;
    G_e(idx) = sqrt( (abs(E_p)^2+abs(E_t)^2)/Z_0 );
end
G_e = 20*log10(G_e);

% Embedded element gain for H-plane scan
G_h = zeros(size(theta));
for idx = 1:length(theta)
    theta_0 = theta(idx);
    phi_0 = 90/180*pi;
    cos_t = cos(theta_0);
    sin_t = sin(theta_0);
    cos_p = cos(phi_0);
    sin_p = sin(phi_0);
    [R_fl,X_fl] = z_fl_dp(theta_0,phi_0,ep_r,a,b,h,l,t);
    nrt_p = sin_p*sinc( k0*(l/2)*sin_t*cos_p/pi ) *...
        sinc( k0*(t/2)*sin_t*sin_p/pi );
    dnrt_p = 1 - 1j*sqrt(ep_r-sin_t^2)*cot(k0*h*sqrt(ep_r-sin_t^2))/cos_t;
    E_p = 1j*I_inc*Z_0*l*(2*R_g/(R_fl+R_g+1j*(X_fl+X_g))) * nrt_p/dnrt_p;
    nrt_t = cos_t*cos_p*sinc( k0*(l/2)*sin_t*cos_p/pi ) *...
        sinc( k0*(t/2)*sin_t*sin_p/pi );
    dnrt_t = 1 - 1j*ep_r*cos_t/sqrt(ep_r-sin_t^2)*cot(k0*h*sqrt(ep_r-sin_t^2));
    E_t = -1j*I_inc*Z_0*l*(2*R_g/(R_fl+R_g+1j*(X_fl+X_g))) * nrt_t/dnrt_t;
    G_h(idx) = sqrt( (abs(E_p)^2+abs(E_t)^2)/Z_0 );
end
G_h = 20*log10(G_h);

% Embedded element gain for D-plane scan
G_d = zeros(size(theta));
for idx = 1:length(theta)
    theta_0 = theta(idx);
    phi_0 = 45/180*pi;
    cos_t = cos(theta_0);
    sin_t = sin(theta_0);
    cos_p = cos(phi_0);
    sin_p = sin(phi_0);
    [R_fl,X_fl] = z_fl_dp(theta_0,phi_0,ep_r,a,b,h,l,t);
    nrt_p = sin_p*sinc( k0*(l/2)*sin_t*cos_p/pi ) *...
        sinc( k0*(t/2)*sin_t*sin_p/pi );
    dnrt_p = 1 - 1j*sqrt(ep_r-sin_t^2)*cot(k0*h*sqrt(ep_r-sin_t^2))/cos_t;
    E_p = 1j*I_inc*Z_0*l*(2*R_g/(R_fl+R_g+1j*(X_fl+X_g))) * nrt_p/dnrt_p;
    nrt_t = cos_t*cos_p*sinc( k0*(l/2)*sin_t*cos_p/pi ) *...
        sinc( k0*(t/2)*sin_t*sin_p/pi );
    dnrt_t = 1 - 1j*ep_r*cos_t/sqrt(ep_r-sin_t^2)*cot(k0*h*sqrt(ep_r-sin_t^2));
    E_t = -1j*I_inc*Z_0*l*(2*R_g/(R_fl+R_g+1j*(X_fl+X_g))) * nrt_t/dnrt_t;
    G_d(idx) = sqrt( (abs(E_p)^2+abs(E_t)^2)/Z_0 );
end
G_d = 20*log10(G_d);

figure(); hold on;
plot(theta/pi*180,G_e,'-k','linewidth',1.1);
plot(theta/pi*180,G_h,'-.r','linewidth',1.1);
plot(theta/pi*180,G_d,'--b','linewidth',1.1);
xlabel('$\theta$ [deg]','interpreter','latex','fontsize',12);
ylabel('Gain [dBi]','interpreter','latex','fontsize',12);
legend('E-plane','H-plane','D-plane','interpreter','latex','fontsize',10.5);
ylim([-25,10]);
grid on;
set(gca,'LooseInset',get(gca,'TightInset'));