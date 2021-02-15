[X,Y] = meshgrid(x,y);
diff = H_fit - H_fit_ZB;

x_hd = linspace(0,0.7,1000);
y_hd = x_hd;
[X_hd,Y_hd] = meshgrid(x_hd,y_hd);
diff_hd = interp2(X,Y,diff,X_hd,Y_hd);
X_hd(X_hd<0.001 | X_hd>0.699) = NaN;
Y_hd(Y_hd<0.001 | Y_hd>0.499) = NaN;
diff_zero_X = X_hd(abs(diff_hd)<1E-5);
diff_zero_Y = Y_hd(abs(diff_hd)<1E-5);

figure(2)
scatter_red = plot(X_hd(diff_hd>=0),Y_hd(diff_hd>=0),'w'); %WZ more stable
hold on
text(0.3,0.15,'Wurtzite','Color','k','FontSize',25)
scatter_blue = plot(X_hd(diff_hd<0),Y_hd(diff_hd<0),'w'); %ZB more stable
text(0.1,0.4,'Zinc Blende','Color','k','FontSize',25)
h = plot(diff_zero_X,diff_zero_Y,'k','LineWidth',2);
x2 = linspace(0,0.7);
y2 = 1-x2;
plot(x2,y2,'k','LineWidth',2)
uistack(h,'top')

shading flat
grid off
box on
hold on
axis equal
view(0,90) %show plot in 2D
xlabel('Indium Mole Fraction')
ylabel('Boron Mole Fraction')
set(gca,'xtick',[0:0.1:0.7])
set(gca,'ytick',[0:0.1:0.5])
% title('Red = WZ more stable, Blue = ZB more stable')
set(gca,'FontSize',15,'Layer','top')
xlim([0 0.7])
ylim([0 0.5])

%Black Cover-Up
[X,Y] = meshgrid(linspace(0,0.7,1000),linspace(0,0.7,1000));
Z2 = ones(size(X)).*0.001;
Z2(X+Y<1) = NaN;
q2 = surf(X,Y,Z2,'FaceColor','k');
set(q2,'edgecolor','none')
