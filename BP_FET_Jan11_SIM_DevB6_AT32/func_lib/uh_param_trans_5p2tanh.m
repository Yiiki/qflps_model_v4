function beta=uh_param_trans_5p2tanh(b)
% b0=beta(4)-beta(3);
% b1=2.*beta(2).^-1;
% b2=2.*beta(2).^-1;
% b3=beta(1);
% b4=beta(3)+beta(4);
% b=[b0,b1,b2,b3,b4];

b4=b(5);b3=b(4);b2=b(3);b1=b(2);b0=b(1);
if b1==b2

    beta4=0.5.*(b0+b4);
    beta3=0.5.*(b4-b0);
    beta2=2.*b1.^-1;
    beta1=b3;
    
    beta=[beta1 beta2 beta3 beta4];
    
else
    error('not matched b1 and b2 !')
end

end