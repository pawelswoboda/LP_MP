function [ ] = Tomo2zimpl(degs,x,file)

n = numel(x);
l = full(max(x))+1;
d = round(sqrt(n));

assert(l > 1);
assert(d^2 == n);

G = gradient_discrete_4(d);
A = tomo_parallel_beam_binary(d,degs,1);
b = A*x;

m = size(A,1);
e = size(G,1);

f=fopen(file,'w');

fprintf(f,'set V:={1..%d};\n',n);
fprintf(f,'set E:={');

  function []=getGrad(k,sep)
    idx = find(abs(G(k,:)));
    assert(numel(idx)==2);
    fprintf(f,'%s<%d,%d>',sep,idx(1),idx(2));
  end

getGrad(1,'')
for i=2:e
  getGrad(i,',');
end

fprintf(f,'};\n');
fprintf(f,'set I:={1..%d};\n',m);
s1 = 'set P[I]:='; s2 = 'param b[I]:=';

  function [] = getProj(k,sep)
    idx = find(A(k,:));
    assert(numel(idx)>0);
    s1 = [s1 sprintf('%s<%d> {%d',sep,k,idx(1))];
    for j=2:numel(idx)
      s1 = [s1 sprintf(',%d',idx(j))];
    end
    s1 = [s1 '}'];
    s2 = [s2 sprintf('%s<%d> %d',sep,k,b(k))];
  end

getProj(1,'');
for i=2:m
  getProj(i,',');
end

fprintf(f,'%s;\n',s1);
fprintf(f,'%s;\n',s2);

fprintf(f,'%s\n','set C:= union<i> in I: P[i];');

fprintf(f,['\n' ...
'var zp[E] real;\n' ...
'var zm[E] real;\n' ...
'var a[E] binary;\n' ...
sprintf('var x[V] integer >= 0 <= %d',l-1) ';\n\n' ...
'minimize obj: sum<i,j> in E: zp[i,j] +  sum<i,j> in E: zm[i,j];\n' ...
'  subto projection:\n' ...
'    forall <i> in I:\n' ...
'      sum<j> in P[i]: x[j] == b[i];\n' ...
'  subto gradient:\n' ...
'    forall <i,j> in E:\n' ...
'      x[i] - x[j] - zp[i,j] + zm[i,j] == 0;\n' ...
'  subto unary:\n' ...
'    forall <i> in V-C:\n' ...
'      x[i] == 0;\n']);

fprintf(f,['' ...
'  subto consistency:\n' ...
'    forall <i,j> in E:\n' ...
'      zp[i,j] <= a[i,j]*2 and zm[i,j] <= (1-a[i,j])*2;\n']);

fclose(f);

cmd = sprintf('zimpl %s',file);
system(cmd);

end

