clear;

Data = importdata('2DHybridExample.dat');
DataSize = size(Data);
numElements = (sqrt(DataSize(1))-1)^2;

h = 1/sqrt(numElements);
[X,Y] = meshgrid(0:h:1,0:h:1);

for i=1:(sqrt(numElements)+1)
    for j=1:(sqrt(numElements)+1)
        Z(i,j) = Data((sqrt(numElements)+1)*(i-1)+j,3);
    end
end

figure;
contourf(X,Y,Z);