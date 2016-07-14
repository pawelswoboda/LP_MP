function [ A ] = tomo_parallel_beam_binary( d,degs,mask )
%% [ A ] = tomo_parallel_beam_binary( d,degs )
% d    -> dimension
% degs -> degrees 0 ... 179
% mask -> 0/1 (if =1 then the matrix only measures pixel inside the circle)
%
%%
  if( nargin < 3 ), mask = 0; end

  n = d^2;
  A = sparse(d*numel(degs),n);

  [I,J] = meshgrid(1:d);J = flipud(J);
  mp = (d+1)/2; % center
  r = d/2;      % radius
  
  M = ( (I-mp).^2+(J-mp).^2 ).^(0.5) <= r;
  M = reshape(M,n,1); 
  M = (M==1);  
  
  for l=1:numel(degs)
    % degree to radian
    deg = degs(l)*(pi/180);
    
    % calculate projection direction p
    % and the orthogonal vector t to p 
    p = [cos(deg);sin(deg)];
    t = [-p(2);p(1)];
    
    % check whether the degree is feasible and
    if( degs(l) <= 180 )
      % n1,n2 are the orthogonal projected end points
      [n1,n2]=proj([mp mp],mp-d,t);
    else
      continue;
    end
    
    % unit step
    step = (n2-n1)/(d-1);
  
    % calculate the right hand sides
    yk = zeros(d,1);
    for k=1:d
      x = n1+(k-1)*step;
      if( sign(x(1)) == 0 ), s=1; else s = -sign(x(1)); end
      yk(k) = s*norm(x,2);
    end

    % create projection
    vars = zeros(1,d^2);
    idx = zeros(d,d);
    for k=1:d-1
      Z1 = abs(yk(k)*ones(d,d) - t(1)*I - t(2)*J);
      Z2 = abs(yk(k+1)*ones(d,d) - t(1)*I - t(2)*J);
      
      M1 = (Z1 <= 0.5 + 1e-6);
      M2 = (Z2 <= 0.5 + 1e-6);
     
      M3 = M1.*M2;
      M1(M3==1)=0;
      M1(idx==1)=1;
      idx = M3;
      
      A((l-1)*d+k,:) = reshape(sparse(M1)',1,d^2);
      vars = vars + A((l-1)*d+k,:);
      if( k == d-1 )
        A((l-1)*d+k+1,:) = reshape(sparse(M2)',1,d^2);
        vars = vars + A((l-1)*d+k+1,:);
      end
    end
    assert(max(vars)==1);    % no pixel is measured twice
    assert(min(vars(M))==1); % every pixel is measured
  end
   
  % set every entry outside the circe to zero 
  if( mask == 1 ) 
    M = ( (I-mp).^2+(J-mp).^2 ).^(0.5) <= r;
    M = reshape(M,n,1); 
    M = (M==0);
    for i=1:size(A,1)
      A(i,M)=0;
    end
  end   

end

function [n1,n2]=proj(mp,r,t)
  mp = (mp*t)*t;
  n1 = mp + r*t;
  n2 = mp - r*t;
end