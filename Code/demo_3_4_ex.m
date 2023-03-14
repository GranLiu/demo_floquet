% -------------------------------------------------------------------------
% Extend Fig.3.4 & Fig.3.5 in Phased Array Antennas.
% 
% Consider printed dipoles with sinc distribution in x and dirac delta in
% y.
%
% Yongxi Liu, Xi'an Jiaotong University, 2023-03.
% -------------------------------------------------------------------------
clc;
clear;
% close all;

addpath("./func");

%% simulation params
f = 2e9;
w = 2*pi*f;
lbd = 3e8/f;
k0 = 2*pi/lbd;
a = 0.6*lbd;
b = 0.6*lbd;
l = 0.4*lbd;
t = 0.1*lbd;
h = 0.05*lbd;
ep_r = 2.0;
ep_0 = 8.854*1e-12;
mu_0 = 1.256*1e-6;
% m = -20:1:20;
m = -100:1:100;
n = m;
d_ang = 0.1;
theta = (d_ang:d_ang:90)/180*pi;

%% E-plane scan
phi = 0/180*pi;

Z_fl_e = zeros(length(theta),1);
for idx_i = 1:length(theta)
    kx0 = k0*sin(theta(idx_i))*cos(phi);
    ky0 = k0*sin(theta(idx_i))*sin(phi);
    kxmn = zeros(length(m),length(n));
    kymn = zeros(length(m),length(n));
    for i = 1:length(m)
        for j = 1:length(n)
            kxmn(i,j) = kx0+2*m(i)*pi/a;
            kymn(i,j) = ky0+2*n(j)*pi/a;
        end
    end

    % if kz is pure imaginary, add minus sign
    kzmn_p = sqrt(k0^2-kxmn.^2-kymn.^2);
    kzmn_m = sqrt(ep_r*k0^2-kxmn.^2-kymn.^2);
    kzmn_p(imag(kzmn_p)~=0) = -kzmn_p(imag(kzmn_p)~=0);
    kzmn_m(imag(kzmn_m)~=0) = -kzmn_m(imag(kzmn_m)~=0);
    
    Y_te_p = kzmn_p/(w*mu_0);
    Y_te_m = kzmn_m/(w*mu_0);
    Y_tm_p = w*ep_0./kzmn_p;
    Y_tm_m = w*ep_0*ep_r./kzmn_m;
    y_te = Y_te_p-1j*Y_te_m.*cot(kzmn_m*h);
    y_tm = Y_tm_p-1j*Y_tm_m.*cot(kzmn_m*h);
%     Z_fl_e(idx_i) = sum( (kymn.^2./y_te+kxmn.^2./y_tm).*...
%         ( sinc(kxmn*l/2/pi) .* sinc(kymn*t/2/pi) ).^2 ./ (k0^2-kzmn_p.^2) ,'all');
    Z_fl_e(idx_i) = sum( (kymn.^2./y_te+kxmn.^2./y_tm) ./ (k0^2-kzmn_p.^2).*...
        ( (cos(kxmn*l/2)-cos(k0*l/2)) ./ (k0^2-kxmn.^2) ).^2 ,'all');
end
% Z_fl_e = l^2/(a*b) * Z_fl_e;
Z_fl_e = 4*k0^2/(a*b) * Z_fl_e;

z_fl_real_e = real(Z_fl_e);
z_fl_imag_e = imag(Z_fl_e);

figure(); hold on;
plot(theta/pi*180,z_fl_real_e);
plot(theta/pi*180,z_fl_imag_e,'--');
ylim([-50,50]);
xlabel("Scan angle (degree)",'interpreter','latex','fontsize',12);
ylabel("$R,X~(\Omega)$",'interpreter','latex','fontsize',12);
legend('R','X','interpreter','latex','fontsize',10);
grid on;
set(gca,'LooseInset',get(gca,'TightInset'));
