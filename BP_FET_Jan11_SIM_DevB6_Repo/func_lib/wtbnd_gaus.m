function out=wtbnd_gaus(faket,pvec,tinq,sig)
% code author: yzy
% date: 2022/12/22
% abstract: (faket, pvec) --> smoothed pvec,
%           at tinq
%           with sig as algorithm parameters 
% algorithm ref: https://doi.org/10.1016/j.chemolab.2003.08.001.

% 1. faket is the time sampling vector with monotonic order;

% 2. pvec is the f sample vector, which should be periodic. Otherwise, some
% patch is suggested, such as linear extrapolation, etc. If an acyclic signal
% is input, the results might be wield at faket(end).

% 3. tinq is the time inquired, which should fulfil min(faket) < tinq < max(faket)

% 4. sig controls the spread of Gaussian function, which tends to a delta
% function if sig -> 0. It is suggested to be 5-10 times of the typical
% spacing of the time sampling vector.


% to extend the available derivatives to n-th order, you need define n-th derivative of
% Gaussian, whose ratio to an Gaussian can be looked up at an ordinary
% mathematical handbook, or online resources such as
% http://www.sci.utah.edu/~gerig/CS7960-S2010/handouts/04%20Gaussian%20derivatives.pdf

if size(pvec,1)==1% unified as collumn vector
pvec=pvec';
end

if size(faket,1)==1% unified as collumn vector
faket=faket';
end

% pvec=[flip(pvec(2:end));pvec;flip(pvec(1:end-1))];
pvec=inflect_extrapolate_func(pvec);
faket=inflect_extrapolate_func(faket);

pvect=pvec-pvec(1);

out=tinq;

for i=1:length(tinq)
    out(i)=inptdat(sig,faket,pvect,tinq(i));
end


out=out+pvec(1);



end

function yinv=inflect_extrapolate_func(y)

yinv=[2.*y(1)-flip(y(2:end));y;2.*y(end)-flip(y(1:end-1))];

end

function y=inptdat(sig,faket,pvec,tinq)
% y=integral(@(t) gaus(sig,tinq-t).*interp1(faket,pvec,t),min(faket),max(faket));
uvec=(faket-tinq).*(sqrt(2).*sig).^-1;
amu=diff(pvec).*diff(faket).^-1;
bmu=pvec(1:end-1)-amu.*uvec(1:end-1).*(sqrt(2).*sig);
erfvec=erf111(uvec);
exp2vc=exp(-uvec.^2);
y=0.5.*bmu'*diff(erfvec)-sig.*(2.*pi).^-0.5.*amu'*diff(exp2vc);

% rectfiying
y=sign(faket(end)-faket(1)).*y;
end

% function y=gaus(sig,t)
% y=sig.^-1.*sqrt(2.*pi).^-1.*exp(-0.5.*sig.^-2.*t.^2);
% end
