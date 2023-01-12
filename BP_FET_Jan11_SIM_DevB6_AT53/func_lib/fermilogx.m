function y=fermilogx(beta,x)
% y=beta(1).*(exp((x-beta(2)).*beta(3).^-1)+1).^-1.*log_exp_plus(-(x-beta(4)).*beta(5).^-1);
y=beta(1).*(exp((x-beta(2)).*beta(3).^-1)+1).^-1.*exp(-x.*beta(4).^-1);
end