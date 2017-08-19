Data = importdata('2DExample.dat');
numElements = 16;

h = 1/sqrt(numElements);
[X,Y] = meshgrid(0:h:1,0:h:1);

for i=1:(sqrt(16)+1)
    for j=1:(sqrt(16)+1)
        Z(i,j) = Data((sqrt(16)+1)*(i-1)+j,3);
    end
end

figure;
surf(X,Y,Z);