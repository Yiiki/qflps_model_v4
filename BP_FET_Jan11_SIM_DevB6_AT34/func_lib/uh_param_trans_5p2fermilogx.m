function beta=uh_param_trans_5p2fermilogx(b)
% b0=0;
% b1=beta(3).^-1;
% b2=-beta(4).^-1;
% b3=beta(2);
% b4=beta(1).*exp(-beta(4).^-1.*beta(2));
% b=[b0,b1,b2,b3,b4];

b4=b(5);b3=b(4);b2=b(3);b1=b(2);b0=b(1);
if b0==0

    beta3=b1.^-1;
    beta4=-b2.^-1;    
    beta2=b3;
    beta1=b4.*exp(-b2.*b3);
    
    beta=[beta1 beta2 beta3 beta4];
    
else
    error('not vanished b0 !')
end
end