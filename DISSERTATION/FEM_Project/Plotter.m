ErrorData = importdata('1DError.dat');
Data = importdata('1DExample.dat');

x = 0:0.01:1;

figure;
plot(Data(:,1), Data(:,2));
hold on
plot(x, (sin(x.*pi)).^2);
legend('FEM','Actual');
xlabel('x');
ylabel('u(x)');

numElements = 1:1:30;
h =  1./numElements;

figure;
loglog(h, h.^(2),'*');

figure;
semilogy(ErrorData(:,1), ErrorData(:,2));
xlabel('$\sqrt[d]{Dofs}$','interpreter', 'latex');
ylabel('Error');