clear all;
Nex = 1280;
data = load(['N',num2str(Nex),'_uv.dat']);
ex_u = data(:,1);
ex_v = data(:,2);
ex_u = reshape(ex_u, Nex+1, Nex+1);
ex_v = reshape(ex_v, Nex+1, Nex+1);

Nmesh = 7;
for i = 1:Nmesh
  n = 10*2^(i-1);
  data = load(['N',num2str(n),'_uv.dat']);
  num_u = data(:,1);
  num_v = data(:,2);
  ex_u_part = ex_u(1:Nex/n:end, 1:Nex/n:end);
  ex_v_part = ex_v(1:Nex/n:end, 1:Nex/n:end);
  ex_u_part = reshape(ex_u_part, (n+1)^2, 1);
  ex_v_part = reshape(ex_v_part, (n+1)^2, 1);
  err(1,i) = norm(num_u-ex_u_part, 1)/n^2;
  err(2,i) = norm(num_v-ex_v_part, 1)/n^2;
  err(3,i) = norm(num_u-ex_u_part, 2)/n;
  err(4,i) = norm(num_v-ex_v_part, 2)/n;
  err(5,i) = norm(num_u-ex_u_part, 'inf');
  err(6,i) = norm(num_v-ex_v_part, 'inf');
end
order = zeros(6, Nmesh-1);
for i = 1:Nmesh-1
  order(:,i) = log2(err(:,i)./err(:,i+1));
end

%0.5897 &  -     & 0.1332 &  -
%0.2852 & 1.0477 & 0.0640 & 1.0568
%0.1386 & 1.0409 & 0.0308 & 1.0574
%0.0669 & 1.0504 & 0.0148 & 1.0585
%0.0317 & 1.0763 & 0.0069 & 1.0992
%0.0134 & 1.2403 & 0.0029 & 1.2287
%0.0046 & 1.5346 & 0.0010 & 1.5527

%function [i,j] = e2ij(N, e)
%j = e/(2*N);
%i = e - 2*N*j;
%end

%function e = ij2e(N, i, j)
%e = 2*N*j + i;
%end

