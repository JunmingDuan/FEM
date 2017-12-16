hehe = load('160.dat');

n = sqrt(size(hehe,1)/2);
hold on;
err = reshape(hehe(1:n*n)-hehe(n*n+1:2*n*n), n, n);
surf(linspace(0,1,n),linspace(0,1,n),err');
shading interp;
colorbar('EastOutside')
axis([0, 1, 0, 1]);
view(0,90);
