function xi=ap_sol_sft_NIXS(X,a,b,xi_n,xi_p,sa,ta,sb,tb)
%     xi=((max(0,xi_p)>xi_n).*a.*xi_n+(min(0,xi_n)<xi_p).*b.*xi_p)...
%         ./(1+(max(0,xi_p)>xi_n).*a+(min(0,xi_n)<xi_p).*b);
    xi = (fq(xi_n,xi_p,sa,ta,b).*a.*xi_n+fq(-xi_p,-xi_n,sb,tb,a).*b.*xi_p)...
        ./(1+fq(xi_n,xi_p,sa,ta,b).*a+fq(-xi_p,-xi_n,sb,tb,a).*b);

    for i=1:X
    % One-Shot Newton-Iteration
    xi = xi - F(a,b,xi_n,xi_p,xi).*dFdx(a,b,xi_n,xi_p,xi).^-1;
    end
end

function z = F(a,b,xi_n,xi_p,xi)
    z = a.*log_exp_plus(xi-xi_n) - b.*log_exp_plus(xi_p-xi) + xi;
end

function z = dFdx(a,b,xi_n,xi_p,xi)
    z = a.*(1+exp(-xi+xi_n)).^-1 + b.*(1+exp(-xi_p+xi)).^-1 + 1;
end

function  z = fq(x,y,s,t,a)
z = ( exp((x - s.*log_exp_plus(a.*(a+1).^-1.*y.*s.^-1)).*t.^-1) + 1 ).^-1;
end

function y=log_exp_plus(x)
% numerical limitation is set by the exp(-inf)
bdn=-6;
y(x<bdn)=log_exp_plus_expand(x(x<bdn));
y(x>=bdn)=log(exp(x(x>=bdn))+1);

    function y=log_exp_plus_expand(x)
        y=exp(x)-exp(2.*x)./2+exp(3.*x)./3-exp(4.*x)./4+...
            exp(5.*x)./5-exp(6.*x)./6+exp(7.*x)./7-exp(8.*x)./8;
    end
end