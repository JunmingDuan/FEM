clear all;
Nex = 640;
data = load(['N',num2str(Nex),'_uv.dat']);
ex_u = data(:,1);
ex_v = data(:,2);
ex_u = reshape(ex_u, Nex+1, Nex+1);
ex_v = reshape(ex_v, Nex+1, Nex+1);

Nmesh = 6;
for i = 2:Nmesh
  n = 10*2^(i-1);
  data = load(['N',num2str(n),'_uv.dat']);
  num_u = data(:,1);
  num_u(1)
  num_v = data(:,2);
  num_u = reshape(num_u, n+1, n+1);
  num_v = reshape(num_v, n+1, n+1);
  ex_u_part = ex_u(1:Nex/n:end, 1:Nex/n:end);
  ex_v_part = ex_v(1:Nex/n:end, 1:Nex/n:end);
  err(1,i) = norm(num_u-ex_u_part, 'fro')/n;
  err(2,i) = norm(num_v-ex_v_part, 'fro')/n;
end
order = zeros(2, 6);
for i = 2:Nmesh-1
  order(:,i) = log2(err(:,i)./err(:,i+1));
end
order

