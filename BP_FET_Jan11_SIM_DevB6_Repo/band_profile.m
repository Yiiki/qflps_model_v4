xtest=linspace(-3,3,1e2);
ytest=erf(-xtest-1);
ytest2=erf(-xtest+1);
figure
tiledlayout(1,1,'TileSpacing','compact')
plot(xtest,ytest2,xtest,ytest,'LineWidth',4)
set(gca,'XColor', 'none','YColor','none')
set(gca, 'color', 'none');