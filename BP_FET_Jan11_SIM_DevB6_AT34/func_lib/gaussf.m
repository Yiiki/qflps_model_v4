function y=gaussf(b,x)
% A, mu, sig, y0
y = b(4)+ b(1).*exp(-((x-b(2))./b(3)).^2);

end