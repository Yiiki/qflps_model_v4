function y=dilog_exp_plus(x)

% numerical limitation is set by the exp(-inf)
bdn=-6;
y(x<bdn)=dilog_exp_plus_expand(x(x<bdn));
y(x>=bdn)=-dilog(1+exp(x(x>=bdn)));
% y(isnan(y))=x(isnan(y));

    function y=dilog_exp_plus_expand(x)
        y=exp(x)-exp(2.*x)./4+exp(3.*x)./9-exp(4.*x)./16+exp(5.*x)./25-exp(6.*x)./36+exp(7.*x)./49-exp(8.*x)./64;
    end
end