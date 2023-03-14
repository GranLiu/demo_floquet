% -------------------------------------------------------------------------
% Realize Fig.3.4 & Fig.3.5 in Phased Array Antennas.
% 
% Yongxi Liu, Xi'an Jiaotong University, 2023-03.
% -------------------------------------------------------------------------
clc;
clear;
close all;

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
% m = -100:1:100;
m = -20:1:20;
n = m;
d_ang = 0.1;
theta = (0.1:d_ang:90)/180*pi;

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
    Z_fl_e(idx_i) = sum( l^2/(a*b)*(kymn.^2./y_te+kxmn.^2./y_tm).*...
        ( sinc(kxmn*l/2/pi) .* sinc(kymn*t/2/pi) ).^2 ./ (k0^2-kzmn_p.^2) ,'all');
end

z_fl_real_e = real(Z_fl_e);
z_fl_imag_e = imag(Z_fl_e);

figure(); hold on;
plot(theta/pi*180,z_fl_real_e);
plot(theta/pi*180,z_fl_imag_e+136,'--');
ylim([-50,50]);
xlabel("Scan angle (degree)",'interpreter','latex','fontsize',12);
ylabel("$R,X~(\Omega)$",'interpreter','latex','fontsize',12);
legend('R','X+136','interpreter','latex','fontsize',10);
grid on;
set(gca,'LooseInset',get(gca,'TightInset'));

%% H-plane scan
phi = 90/180*pi;
Z_fl_h = zeros(length(theta),1);
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
    Z_fl_h(idx_i) = sum( l^2/(a*b)*(kymn.^2./y_te+kxmn.^2./y_tm).*...
        ( sinc(kxmn*l/2/pi) .* sinc(kymn*t/2/pi) ).^2 ./ (k0^2-kzmn_p.^2) ,'all');
end

z_fl_real_h = real(Z_fl_h);
z_fl_imag_h = imag(Z_fl_h);

figure(); hold on;
plot(theta/pi*180,z_fl_real_h);
plot(theta/pi*180,z_fl_imag_h+136,'--');
ylim([0,20]);
xlabel("Scan angle (degree)",'interpreter','latex','fontsize',12);
ylabel("$R,X~(\Omega)$",'interpreter','latex','fontsize',12);
legend('R','X+136','interpreter','latex','fontsize',10);
grid on;
set(gca,'LooseInset',get(gca,'TightInset'));

%% D-plane scan
phi = 45/180*pi;
Z_fl_d = zeros(length(theta),1);
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
    Z_fl_d(idx_i) = sum( l^2/(a*b)*(kymn.^2./y_te+kxmn.^2./y_tm).*...
        ( sinc(kxmn*l/2/pi) .* sinc(kymn*t/2/pi) ).^2 ./ (k0^2-kzmn_p.^2) ,'all');
end

z_fl_real_d = real(Z_fl_d);
z_fl_imag_d = imag(Z_fl_d);

figure(); hold on;
plot(theta/pi*180,z_fl_real_d);
plot(theta/pi*180,z_fl_imag_d+136,'--');
ylim([0,20]);
xlabel("Scan angle (degree)",'interpreter','latex','fontsize',12);
ylabel("$R,X~(\Omega)$",'interpreter','latex','fontsize',12);
legend('R','X+136','interpreter','latex','fontsize',10);
grid on;
set(gca,'LooseInset',get(gca,'TightInset'));

%% reflection coefficient
% power wave model
z_0 = conj(Z_fl_e(1));

ref_co_e = abs( (Z_fl_e-conj(z_0))./(Z_fl_e+z_0) );
ref_co_h = abs( (Z_fl_h-conj(z_0))./(Z_fl_h+z_0) );
ref_co_d = abs( (Z_fl_d-conj(z_0))./(Z_fl_d+z_0) );

% tranmission line model
% z_0 = Z_fl_e(1);
% 
% ref_co_e = abs( (Z_fl_e-z_0)./(Z_fl_e+z_0) );
% ref_co_h = abs( (Z_fl_h-z_0)./(Z_fl_h+z_0) );
% ref_co_d = abs( (Z_fl_d-z_0)./(Z_fl_d+z_0) );

figure(); hold on;
plot(theta/pi*180,ref_co_e);
plot(theta/pi*180,ref_co_h,'--');
plot(theta/pi*180,ref_co_d,'-o','MarkerSize',1);
ylim([0,1]);
xlabel("Scan angle (degree)",'interpreter','latex','fontsize',12);
ylabel("Reflection coefficienct",'interpreter','latex','fontsize',12);
legend('$\phi=0^{\circ}$','$\phi=45^{\circ}$','$\phi=90^{\circ}$','interpreter','latex','fontsize',10);
grid on;
set(gca,'LooseInset',get(gca,'TightInset'));
