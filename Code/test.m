figure(); hold on;
plot(theta/pi*180,z_fl_real);
plot(theta/pi*180,z_fl_imag+136,'--');
% ylim([-50,50]);
title('$m=-10:1:10$','interpreter','latex','fontsize',12);
xlabel("Scan angle (degree)",'interpreter','latex','fontsize',12);
ylabel("$R,X~(\Omega)$",'interpreter','latex','fontsize',12);
legend('R','X','interpreter','latex','fontsize',10);
grid on;
set(gca,'LooseInset',get(gca,'TightInset'));

x = (-20:0.001:20);
y = cot(1j * x);

figure();
plot(x,real(y));

figure();
plot(x,imag(y));