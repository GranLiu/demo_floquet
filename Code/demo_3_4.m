% -------------------------------------------------------------------------
% Realize Fig. 3.4 in Phased Array Antennas.
% 
% 
% Yongxi Liu, Xi'an Jiaotong University, 2023-03.
% -------------------------------------------------------------------------
clc;
clear;
% close all;

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
m = -10:1:10;
n = m;

%% E-plane scan
phi = 0/180*pi;
theta = (0:0.5:90)/180*pi;
Z_fl = zeros(length(theta),1);
for idx_i = 1:length(theta)
    kx0 = k0*sin(theta(idx_i))*cos(phi);
    ky0 = k0*sin(theta(idx_i))*sin(phi);
    for idx_m = 1:length(m)
        for idx_n = 1:length(n)
            kxmn = kx0+2*m(idx_m)*pi/a;
            kymn = ky0+2*n(idx_n)*pi/a;
            kzmn_p = sqrt(k0^2-kxmn^2-kymn^2);
            kzmn_m = sqrt(ep_r*k0^2-kxmn^2-kymn^2);
            Y_te_p = kzmn_p/(w*mu_0);
            Y_te_m = kzmn_m/(w*mu_0);
            Y_tm_p = w*ep_0/kzmn_p;
            Y_tm_m = w*ep_0*ep_r/kzmn_m;
%             y_te = Y_te_p-1j*Y_te_m*cot(kzmn_m*h);
%             y_tm = Y_tm_p-1j*Y_tm_m*cot(kzmn_m*h);
            y_te = Y_te_p+1j*Y_te_m*cot(kzmn_m*h);
            y_tm = Y_tm_p+1j*Y_tm_m*cot(kzmn_m*h);
            Z_fl(idx_i) = Z_fl(idx_i) + l^2/(a*b)*(kymn^2/y_te+kxmn^2/y_tm)*...
                ( sinc(kxmn*l/2/pi) * sinc(kymn*t/2/pi) )^2/(k0^2-kzmn_p^2);
        end
    end
end

z_fl_real = real(Z_fl);
z_fl_imag = imag(Z_fl);

figure(); hold on;
plot(theta/pi*180,z_fl_real);
plot(theta/pi*180,z_fl_imag,'--');
% ylim([-50,50]);
xlabel("Scan angle (degree)",'interpreter','latex','fontsize',12);
ylabel("$R,X~(\Omega)$",'interpreter','latex','fontsize',12);
legend('R','X','interpreter','latex','fontsize',10);
grid on;
set(gca,'LooseInset',get(gca,'TightInset'));
