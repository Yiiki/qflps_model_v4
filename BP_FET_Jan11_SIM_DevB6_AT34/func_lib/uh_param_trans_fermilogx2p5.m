function b=uh_param_trans_fermilogx2p5(beta)
b0=0;
b1=beta(3).^-1;
b2=-beta(4).^-1;
b3=beta(2);
b4=beta(1).*exp(-beta(4).^-1.*beta(2));
b=[b0,b1,b2,b3,b4];
end