clearvars -except H_fit X0 Y0
data = xlsread('../data/ZB_BInGaN_data_KG.xlsx','Master List - Full','C3:E89'); % Zinc Blende
B_y = data(:,1);
In_x = data(:,2);
a = data(:,3);

% Lattice Parameters [A] (Zinc Blende)
a_BN = 3.621108479;         % From Logan (for our functional)
a_InN = 5.014327652;        % From Logan (for our functional)
a_GaN = 4.51879135581962;   % From Logan (for our functional)

% Get fit parameters (same results as gnuplot)
tbl = table(In_x,B_y,a);
modelfunc = @(p,x) x(:,2).*a_BN + x(:,1).*a_InN + (1-x(:,1)-x(:,2)).*a_GaN + x(:,2).*x(:,1).*p(1) + x(:,1).*(1-x(:,1)-x(:,2)).*p(2) + x(:,2).*(1-x(:,1)-x(:,2)).*p(3); %bowing model
p0 = [1 1 1];
model = fitnlm(tbl,modelfunc,p0);
p = model.Coefficients.Estimate
p_standardError = model.Coefficients.SE;
p_percentError = abs(p_standardError./p).*100

% Plot surface of best fit and horizontal plane at a_GaN
figure(4)
[X,Y] = meshgrid(linspace(0,0.7,2000),linspace(0,0.6,2000));
a_fit3 = Y.*a_BN + X.*a_InN + (1-X-Y).*a_GaN + Y.*X.*p(1) + X.*(1-X-Y).*p(2) + Y.*(1-X-Y).*p(3);

%Relative Error Plot (with color bar)
figure(4);
relErr = (a_fit3-a_GaN)/a_GaN;
relErr(Y>1-X) = NaN;
surf(X,Y,relErr);
hold on
[~,h] = contour(X,Y,relErr,-0.10:0.01:0.10,'--k');
h.ContourZLevel = max(max(relErr))+0.002;
xlabel('Indium Mole Fraction')
ylabel('Boron Mole Fraction')
zlabel('(a - a_{GaN})/a_{GaN}')
title('Lattice Parameter (a): Relative Difference with GaN')
shading flat
grid off
box on
c = colorbar;
set(c,'FontSize',15)
hold on
view(0,90) %show plot in 2D
set(gca,'FontSize',15,'Layer','top')
colorcet('D1')
caxis([-0.05 0.05])
set(c,'YTick',[-0.05:0.01:0.05],'TickLabels',{'-5%', '-4%', '-3%', '-2%', '-1%', ' 0%', ' 1%', ' 2%', ' 3%', ' 4%', ' 5%'},'TickLength',0)

%Plot black line along where relative error is zero
relErr0 = abs(relErr) < 4E-7; %4E-7 for ZB
X0_ZB = X(relErr0);
Y0_ZB = Y(relErr0);
plot3(X0_ZB,Y0_ZB,zeros(size(X0_ZB)),'k')
axis equal
xlim([0 0.7]) % 0.7
ylim([0 0.5]) % 0.5
set(findobj(gca, 'Type', 'Line', 'Linestyle', '-'), 'LineWidth', 2); %Make line thicker

%White Cover-Up (newer version)
[X,Y] = meshgrid(linspace(0,0.7,1000),linspace(0,0.7,1000));
Z = ones(size(X)).*max(max(relErr))+0.001;
Z(Y<X & Y>0.5*X & Y<1-X) = NaN;
q = surf(X,Y,Z,'FaceColor','w');
set(q,'edgecolor','none')
