function plot2d(n)
P = load(['P',num2str(n),'.dat']);
TRI = load(['TRI',num2str(n),'.dat']);

data = load(['N',num2str(n),'_uv.dat']);
N = sqrt(size(data,1));
x = linspace(0,10,N);
y = linspace(0,10,N);

num_u = data(:,1);
num_v = data(:,2);
exa_u = data(:,3);
exa_v = data(:,4);
norm(num_u-exa_u)/N
norm(num_v-exa_v)/N

figure(1)
surf(x, y, reshape(num_u, N, N)', 'EdgeColor', 'interp');
colormap jet
colorbar
caxis([0,0.2]);
xlabel('x-displament');
%view(0,90);

figure(2)
surf(x, y, reshape(num_v, N, N)', 'EdgeColor', 'interp');
colormap jet
colorbar
xlabel('y-displament');

data = load(['N',num2str(n),'_es.dat']);
num_ex = data(:,1);
num_ey = data(:,2);
num_exy = data(:,3);
num_sx = data(:,4);
num_sy = data(:,5);
num_sxy = data(:,6);

%norm(num_ex-exa_ex)/N
%norm(num_ey-exa_ey)/N
%norm(num_exy-exa_exy)/N
%norm(num_sx-exa_sx)/N
%norm(num_sy-exa_sy)/N
%norm(num_sxy-exa_sxy)/N

%norm((num_ex-exa_ex)./exa_ex)/N
%norm((num_ey-exa_ey)./exa_ey)/N
%norm((num_exy-exa_exy)./exa_exy)/N
%norm((num_sx-exa_sx)./exa_sx)/N
%norm((num_sy-exa_sy)./exa_sy)/N
%norm((num_sxy-exa_sxy)./exa_sxy)/N

%figure(3)
%patch('Faces', TRI80, 'Vertices', P40, 'FaceVertexCData', exa_ex, 'FaceColor', 'interp', 'EdgeColor', 'interp');
%patch('Faces', TRI, 'Vertices', P, 'FaceVertexCData', num_ex, 'FaceColor', 'flat', ...
%'EdgeColor', 'none');
%axis([0,0.1,0,0.1]);
%colorbar;
%axis equal;
%xlabel('strain x');

end

