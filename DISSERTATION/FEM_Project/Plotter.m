clear;

Error1 = importdata('1DError_p1.dat');
Error2 = importdata('1DError_p2.dat');
Error3 = importdata('1DError_p3.dat');
Error4 = importdata('1DError_p4.dat');
Data = importdata('1DExample.dat');

x = 0:0.01:1;

figure;
plot(x, (sin(x.*pi)).^2);
hold on
plot(Data(:,1), Data(:,2));
legend('FEM','Actual');
xlabel('x');
ylabel('u(x)');
title('Plot of $u(x) = \sin(\pi x)^2$ for $h=0.0625$, $p=1$','interpreter','latex');

figure;
loglog(Error1(:,1), Error1(:,2));
hold on;
loglog(Error2(:,1), Error2(:,2));
loglog(Error3(:,1), Error3(:,2));
loglog(Error4(:,1), Error4(:,2));
xline = 2e-3;
line([32 64],[xline xline/4],'Color','k');
line([32 64],[xline xline],'Color','k');
line([64 64],[xline xline/4],'Color','k');
text(32+13, xline*2,'1');
text(64+3, xline/2,'2');
xline = 2.5e-5;
line([64 128],[xline xline/8],'Color','k');
line([64 128],[xline xline],'Color','k');
line([128 128],[xline xline/8],'Color','k');
text(64+26, xline*2,'1');
text(128+6, xline/3,'3');
xline = 3e-7;
line([96 192],[xline xline/16],'Color','k');
line([96 192],[xline xline],'Color','k');
line([192 192],[xline xline/16],'Color','k');
text(96+39, xline*2,'1');
text(192+9, xline/4,'4');
xline = 3e-9;
line([128 256],[xline xline/32],'Color','k');
line([128 256],[xline xline],'Color','k');
line([256 256],[xline xline/32],'Color','k');
text(128+52, xline*2,'1');
text(256+12, xline/5,'5');

xlabel('$\sqrt[d]{Dofs}$','interpreter','latex');
ylabel('Error');
legend('p=1','p=2','p=3','p=4');
title('');