function y=erf111(x)
y=x;

xa=0.46875;
xb=4.0;

xan=-xa;
xbn=-xb;

y(abs(x)<=xa)=AR11(x(abs(x)<=xa));
y(x>=xa&x<xb)=1-BR11(x(x>=xa&x<xb));
y(x>=xb)=1-CR11(x(x>=xb));

y(x<xan&x>xbn)=BR11(-x(x<xan&x>xbn))-1;
y(x<=xbn)=CR11(-x(x<=xbn))-1;

function y=AR11(x)
% p 
p0=3.6767877e00;
p1=-9.7970465e-02;
% q
q0=3.2584593e00;
q1=1.0e00;
% numerator
nu=p0.*x.^0+p1.*x.^2;
% denominator
de=q0.*x.^0+q1.*x.^2;
% y
y=x.*nu.*de.^-1;
end

function y=BR11(x)
% p 
p0=7.3033e-01;
p1=-2.3877e-02;
% q
q0=6.6211e-01;
q1=1.0e00;
% numerator
nu=p0.*x.^0+p1.*x.^1;
% denominator
de=q0.*x.^0+q1.*x.^1;
% y
y=exp(-x.^2).*nu.*de.^-1;
end

function y=CR11(x)
% p 
p0=-1.24368544e-01;
p1=-9.68210364e-02;
% q
q0=4.40917061e-01;
q1=1.0e00;
% numerator
nu=p0.*x.^-0+p1.*x.^-2;
% denominator
de=q0.*x.^-0+q1.*x.^-2;
% y
y=exp(-x.^2).*x.^-1.*(pi.^-0.5+x.^-2.*nu.*de.^-1);
end

end
% function y=AR44(x)
% p0=3.209377589138469472562;
% p1=3.209377589138469472562;
