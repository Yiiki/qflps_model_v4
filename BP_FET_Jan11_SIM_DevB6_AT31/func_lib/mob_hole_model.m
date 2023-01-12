function y=mob_hole_model(b,x)
% y = b(1).*exp(-(x-b(2)).^2.*b(3).^-2).*(1 - b(6).*exp(-(x-b(4)).^2.*b(5).^-2));
y=b(1).*tanh((x-b(2))./b(3))+b(4);
end