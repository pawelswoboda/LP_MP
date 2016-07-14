function [  ] = TomoConverter( degs,x,file )

[~,~,ext] = fileparts(file);
strrep(file,ext,'');

Tomo2zimpl(degs,x,[file '.zpl']);
Tomo2MP(degs,x,[file '_l1.txt']);
Tomo2MPf(degs,x,[file '_feas.txt']);

end

