function [Eopt, Dopt, Vopt, Iopt, Sopt, wopt, demopt, sopt] = parse_CVX_output(X,nvar,feeder,nodes,lines,configs,loads,caps,controllers)

nnode = nodes.nnode;
nline = lines.nline;

Xa = X(1:nvar);
Xb = X(nvar+1:2*nvar);
Xc = X(2*nvar+1:3*nvar);

XX = [Xa Xb Xc]';

Eopt = XX(:,1:nnode);
XX(:,1:nnode) = [];

Dopt = XX(:,1:nnode);
XX(:,1:nnode) = [];

Vopt = sqrt(Eopt).*exp(j*pi/180*Dopt).*nodes.PH;

for k1 = 1:lines.nline
    Iopt(:,k1) = lines.FYpu(:,:,k1)*(Vopt(:,lines.TXnum(k1)) - Vopt(:,lines.RXnum(k1)));
end
Iopt = Iopt.*lines.PH;

Popt = XX(:,1:nline);
XX(:,1:nline) = [];

Qopt = XX(:,1:nline);
XX(:,1:nline) = [];

Sopt = Popt + 1j*Qopt;

uopt = XX(:,1:nnode);
XX(:,1:nnode) = [];

vopt = XX(:,1:nnode);
XX(:,1:nnode) = [];

wopt = uopt + 1j*vopt;

demopt = loads.spu.*(loads.aPQ + loads.aZ.*Eopt).*nodes.PH;

for k1 = 2:nnode
    sopt(:,k1) = sum(Sopt(:,nodes.inmat((nodes.inmat(:,k1) ~= 0),k1)),2) - ...
        sum(Sopt(:,nodes.outmat((nodes.outmat(:,k1) ~= 0),k1)),2);
end

end