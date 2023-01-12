function y=uh_model_5p(b,x)
b4=b(5);b3=b(4);b2=b(3);b1=b(2);b0=b(1);
y=(b4.*exp(b2.*(x-b3))+b0).*(exp(b1.*(x-b3))+1).^-1;
end