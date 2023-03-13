function [R,X] = z_fl_dp(theta,phi,ep_r,a,b,h,l,t)
% -------------------------------------------------------------------------
% Compute the Floquet impedance of an infinite array composed of printed
% dipoles, given the scan angle.
%
% Input
%	theta: scan polar angle, in radian
%	phi: azimuth angle, in radian
%   ep_r: relative permittivity
%   a: length of the cell (period in x dimension), in wavelength
%   b: width of the cell (period in y dimension), in wavelength
%   h: height of the dielectirc substrate, in wavelength
%   l: current length (in x dimension), in wavelength
%   t: current width (in y dimension), in wavelength
% Output
%   R: real part of the Z_fl
%   X: imag part of the Z_fl
%
% Reference
% [1] Phased Array Antennas: Floquet Analysis, Synthesis, BFNs and Active
%     Array Systems, ARUN K. BHATTACHARYYA.
%
% Yongxi Liu, Xi'an Jiaotong University, 2023-03.
% -------------------------------------------------------------------------
%% check inputs
if nargin~= 8
    error("There should be 8 inputs in this func.");
end
numchk = {'numeric'};
poschk = {'positive'};
validateattributes(theta,numchk,{'>=',0,'<=',2*pi})
validateattributes(phi,numchk,{'>=',0,'<=',2*pi})
validateattributes(ep_r,numchk,poschk)
validateattributes(a,numchk,poschk)
validateattributes(b,numchk,poschk)
validateattributes(h,numchk,poschk)
validateattributes(l,numchk,poschk)
validateattributes(t,numchk,poschk)

%% params
ep_0 = 8.854*1e-12;
mu_0 = 1.256*1e-6;
m = -20:1:20;
n = m;
k0 = 2*pi;
w = 3e8*k0;

%% simulation
kx0 = k0*sin(theta)*cos(phi);
ky0 = k0*sin(theta)*sin(phi);
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
Z_fl = sum( l^2/(a*b)*(kymn.^2./y_te+kxmn.^2./y_tm).*...
    ( sinc(kxmn*l/2/pi) .* sinc(kymn*t/2/pi) ).^2 ./ (k0^2-kzmn_p.^2) ,'all');

R = real(Z_fl);
X = imag(Z_fl);
end