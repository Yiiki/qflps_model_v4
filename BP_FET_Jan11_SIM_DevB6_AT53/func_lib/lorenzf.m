function y=lorenzf(b,x)
y=b(4)+b(1).*(1 + ((x-b(2)).*b(3).^-1).^2).^-1;
end