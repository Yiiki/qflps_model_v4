function y=phi_model(b,x)

% % parameters
% a=phi_parm(1);
% b=phi_parm(2);
% c=phi_parm(3);
% d=phi_parm(4);
% 
% % formula
% y=c.*tanh((x-a).*b.^-1)+d;

% formula
% y=b(1).*exp(a.*log(log(1+exp((x-b(2)).*b(3).^-1)))) + b(4);
y=b(3).*tanh((x-b(1)).*b(2).^-1)+b(4);

end
