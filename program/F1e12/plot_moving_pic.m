function plot_moving_pic(n, f, i)
P = load(['P',num2str(n),'.dat']);
TRI = load(['TRI',num2str(n),'.dat']);

data = load(['N',num2str(n),'_uv_f',num2str(f/1e11,'%d'),'e11.dat']);
N = sqrt(size(data,1));
x = linspace(0,10,N);
y = linspace(0,10,N);

num_u = data(:,1);
num_v = data(:,2);

%clf;
%surf(x, y, reshape(num_u, N, N)', 'EdgeColor', 'interp');
%colormap jet
%colorbar
%xlabel('x-displament');
%caxis([0,0.025]);
%view(0,90);
%print(['N',num2str(n),'_u_f',num2str(i),'.eps'], '-depsc');

%clf;
%surf(x, y, reshape(num_v, N, N)', 'EdgeColor', 'interp');
%colormap jet
%colorbar
%xlabel('y-displament');
%caxis([-15e-3,1e-3]);
%view(0,90);
%print(['N',num2str(n),'_v_f',num2str(i),'.eps'], '-depsc');
%
%data = load(['N',num2str(n),'_es.dat']);
data = load(['N',num2str(n),'_es_f',num2str(f/1e11,'%d'),'e11.dat']);
num_ex = data(:,1);
num_ey = data(:,2);
num_exy = data(:,3);
num_sx = data(:,4);
num_sy = data(:,5);
num_sxy = data(:,6);
%
%clf;
patch('Faces', TRI, 'Vertices', P, 'FaceVertexCData', num_ex, 'FaceColor', 'flat', ...
'EdgeColor', 'none');
colormap jet
colorbar;
axis equal;
xlabel('strain x');
view(0,90);
axis([0,0.1,0,0.1]);
caxis([-6e3,1e4]);
print(['N',num2str(n),'_ex_f',num2str(i),'.eps'], '-depsc');
max(num_ex)
min(num_ex)
%
%%clf;
%patch('Faces', TRI, 'Vertices', P, 'FaceVertexCData', num_ey, 'FaceColor', 'flat', ...
%'EdgeColor', 'none');
%colormap jet
%colorbar;
%axis equal;
%xlabel('strain y');
%view(0,90);
%axis([0,0.1,0,0.1]);
%caxis([-6e3,8e3]);
%print(['N',num2str(n),'_ey_f',num2str(i),'.eps'], '-depsc');
%max(num_ey)
%min(num_ey)
%
%%clf;
%patch('Faces', TRI, 'Vertices', P, 'FaceVertexCData', num_exy, 'FaceColor', 'flat', ...
%'EdgeColor', 'none');
%colormap jet
%colorbar;
%axis equal;
%xlabel('strain xy');
%view(0,90);
%axis([0,0.1,0,0.1]);
%print(['N',num2str(n),'_exy.eps'], '-depsc');
%
%%clf;
%patch('Faces', TRI, 'Vertices', P, 'FaceVertexCData', num_sx, 'FaceColor', 'flat', ...
%'EdgeColor', 'none');
%colormap jet
%colorbar;
%axis equal;
%xlabel('stress x');
%view(0,90);
%axis([0,0.1,0,0.1]);
%print(['N',num2str(n),'_sx.eps'], '-depsc');
%
%%clf;
%patch('Faces', TRI, 'Vertices', P, 'FaceVertexCData', num_sy, 'FaceColor', 'flat', ...
%'EdgeColor', 'none');
%colormap jet
%colorbar;
%axis equal;
%xlabel('stress y');
%view(0,90);
%axis([0,0.1,0,0.1]);
%print(['N',num2str(n),'_sy.eps'], '-depsc');
%
%%clf;
%patch('Faces', TRI, 'Vertices', P, 'FaceVertexCData', num_sxy, 'FaceColor', 'flat', ...
%'EdgeColor', 'none');
%colormap jet
%colorbar;
%axis equal;
%xlabel('stress xy');
%view(0,90);
%axis([0,0.1,0,0.1]);
%print(['N',num2str(n),'_sxy.eps'], '-depsc');
%
end

