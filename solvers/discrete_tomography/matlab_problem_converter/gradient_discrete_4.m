function [ B ] = gradient_discrete_4( d )

op = speye(d)-triu(ones(d,d),1)+triu(ones(d,d),2);
op = op(1:d-1,:);
B = [ kron(speye(d),op) ;
      kron(op,speye(d)) ];

end

