function [] = Tomo2MP(degs,x,file)

n = numel(x);
l = full(max(x))+1;
d = round(sqrt(n));

assert(l > 1);
assert(d^2 == n);

G = gradient_discrete_4(d);
A = tomo_parallel_beam_binary(d,degs,1);
b = A*x;

check = sum(A,1);

f = fopen(file,'w');

fprintf(f,'MARKOV\n');
fprintf(f,'%d\n',n);
fprintf(f,'%d ',l*ones(n,1));
fprintf(f,'\n');

% Pairwise Factors
fprintf(f,'%d\n',size(G,1)+n);
for i=1:size(G,1)
  idx = find(G(i,:));
  fprintf(f,'%d %d %d\n',2,idx(1)-1,idx(2)-1);
end

% Unary Factors
for i=1:n, fprintf(f,'%d %d\n',1,i-1); end
fprintf(f,'\n');

% Pairwise Factor Tabels
[I,J]=meshgrid(0:l-1);
reg = abs(I-J); reg = reg(:);
for i=1:size(G,1)
  fprintf(f,'%d\n',numel(reg));
  fprintf(f,'%.4f ',reg);
  fprintf(f,'\n');
end

% Unary Factor Tables
un = zeros(1,l);
unI = un; unI(2) = Inf; unI(3) = Inf;
for i=1:n
  fprintf(f,'%d\n',l);
  if( check(i) > 0 )
    fprintf(f,'%.4f ',un);
  else
    fprintf(f,'%.4f ',unI);
  end
  fprintf(f,'\n');
end

fprintf(f,'\n');
fprintf(f,'PROJECTIONS\n');
for i=1:size(A,1)
  idx = find(A(i,:));
  assert(numel(idx)>0);
  fprintf(f,'%d ',idx(1)-1);
  for j=2:numel(idx)
    fprintf(f,'+ %d ',idx(j)-1);
  end
  fprintf(f,'= (');
  for j=0:b(i)-1
    fprintf(f,'Inf,');
  end
  fprintf(f,'0)\n');
end


fclose(f);

end

