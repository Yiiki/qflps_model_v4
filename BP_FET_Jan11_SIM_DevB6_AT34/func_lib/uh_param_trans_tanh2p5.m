function b=uh_param_trans_tanh2p5(beta)
b0=beta(4)-beta(3);
b1=2.*beta(2).^-1;
b2=2.*beta(2).^-1;
b3=beta(1);
b4=beta(3)+beta(4);
b=[b0,b1,b2,b3,b4];
end