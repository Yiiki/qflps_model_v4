function ymat=func_elec(phi,phi0,kbv,vgslis,vdslis)

cols=length(vgslis);
rows=length(vdslis);

ymat=zeros(rows,cols);

for j=1:cols
    Vgs = vgslis(j);
    y = arrayfun(@(x) Func_elec(phi,phi0,kbv,Vgs,x),vdslis);
    if size(vdslis,2)~=1
        y=y';
    end
    ymat(:,j)=y;
end

function out=Func_elec(phi,phi0,kbv,Vgs,Vds0)
% it makes sense for Vds > 0
% requires Vgs > Vthn

An=kbv(1).^-1;
Rsd=kbv(2);
Vthn=kbv(3);

Vds = min(Vds0,Vgs-Vthn);

out = An.*(Vgs-Vthn-0.5.*Vds).*(1+Rsd.*An.*(Vgs-Vthn)).^-1.*Vds.*(1+exp((phi-Vds).*phi0.^-1)).^-1;

end
end