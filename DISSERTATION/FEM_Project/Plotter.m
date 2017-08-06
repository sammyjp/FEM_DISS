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

numElements = 6:1:30;
h = 1./numElements;

figure;
loglog(Error1(:,1), Error1(:,2));
hold on;
loglog(Error2(:,1), Error2(:,2));
loglog(Error3(:,1), Error3(:,2));
loglog(Error4(:,1), Error4(:,2));
loglog(numElements,h.^2);
loglog(numElements,h.^3);
loglog(numElements,h.^4);
loglog(numElements,h.^5);
xlabel('$\sqrt[d]{Dofs}$','interpreter','latex');
ylabel('Error');
legend('p=1','p=2','p=3','p=4');
title('');