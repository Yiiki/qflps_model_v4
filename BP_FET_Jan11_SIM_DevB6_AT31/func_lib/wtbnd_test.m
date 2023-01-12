xtest=linspace(-3,3,1e2);
ytest=tanh(xtest);

xctrl=linspace(-3,3,7);
yctrl=tanh(xctrl);
ymaki=makima(xctrl,yctrl,xtest);

sig=1.0e-1;
yrcos=wtbnd_gaus(xctrl,yctrl,xtest,sig);

figure
plot(xtest,ytest,'-',xctrl,yctrl,'o',xtest,yrcos,'--')

figure
plot(xtest,ytest,'-',xctrl,yctrl,'o',xtest,ymaki,'--')

xtest=linspace(-3,3,1e2);
ytest=exp(-xtest.^2);

xctrl=linspace(-3,3,7);
yctrl=exp(-xctrl.^2);
ymaki=makima(xctrl,yctrl,xtest);

sig=1.0e-1;
yrcos=wtbnd_gaus(xctrl,yctrl,xtest,sig);

figure
plot(xtest,ytest,'-',xctrl,yctrl,'o',xtest,yrcos,'--')

figure
plot(xtest,ytest,'-',xctrl,yctrl,'o',xtest,ymaki,'--')