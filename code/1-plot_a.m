clear
data = xlsread('../data/WZ_BInGaN_data_KG.xlsx','Master List - Full','C3:E115'); % Wurtzite
B_y = data(:,1);
In_x = data(:,2);
a = data(:,3);

% Lattice Parameters [A] (Wurtzite)
a_BN = 2.55159;
a_InN = 3.55463;
a_GaN = 3.192667772;

% Get fit parameters (same results as gnuplot)
tbl = table(In_x,B_y,a);
modelfunc = @(p,x) x(:,2).*a_BN + x(:,1).*a_InN + (1-x(:,1)-x(:,2)).*a_GaN + x(:,2).*x(:,1).*p(1) + x(:,1).*(1-x(:,1)-x(:,2)).*p(2) + x(:,2).*(1-x(:,1)-x(:,2)).*p(3); %bowing model
p0 = [1 1 1];
model = fitnlm(tbl,modelfunc,p0);
p = model.Coefficients.Estimate
p_standardError = model.Coefficients.SE;
p_percentError = abs(p_standardError./p).*100

% %Plot surface of best fit and horizontal plane at a_GaN
figure(31)
[X,Y] = meshgrid(linspace(0,0.7,2000),linspace(0,0.6,2000));
a_fit3 = Y.*a_BN + X.*a_InN + (1-X-Y).*a_GaN + Y.*X.*p(1) + X.*(1-X-Y).*p(2) + Y.*(1-X-Y).*p(3);

%Relative Error Plot (with color bar)
figure(31);
relErr = (a_fit3-a_GaN)/a_GaN;
surf(X,Y,relErr);
hold on
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

% Isobars
spacing = 0.01;
for iso = -0.10:spacing:0.10 %For all isobars, use -0.1:spacing:0.1....-0.105 to 0.095....-0.10 to 0.10
    X_iso = X(abs(relErr - iso) < 1E-5); %5E-7 for ZB
    Y_iso = Y(abs(relErr - iso) < 1E-5); %5E-7 for ZB
    plot3(X_iso,Y_iso,ones(size(X_iso)).*max(max(relErr))+0.002,'--k')
end

%Plot black line along where relative error is zero
relErr0 = abs(relErr) < 5E-5; %4E-7 for ZB
X0 = X(relErr0);
Y0 = Y(relErr0);
plot3(X0,Y0,zeros(size(X0)),'k')

%Border around region with data
x = linspace(0,2/3);
axis equal
xlim([0 0.7]) % 0.7
ylim([0 0.5]) % 0.5
set(findobj(gca, 'Type', 'Line', 'Linestyle', '-'), 'LineWidth', 2); %Make line thicker

%White Cover-Up
[X,Y] = meshgrid(linspace(0,0.7,1000),linspace(0,0.7,1000));
Z = ones(size(X)).*max(max(relErr))+0.001;
Z(Y<X & Y>0.5*X) = NaN;
q = surf(X,Y,Z,'FaceColor','w');
set(q,'edgecolor','none')
Z2 = ones(size(X)).*max(max(relErr))+0.003;
Z2(X+Y<1) = NaN;
q2 = surf(X,Y,Z2,'FaceColor','w');
set(q2,'edgecolor','none')

%%
% % Getting data along lattice-matched line and printing to XLSX file
% p_poly = polyfit(X0,Y0,3);
% x_poly = linspace(0,0.7,71);
% y_poly = polyval(p_poly,x_poly);
% T = table([x_poly]',[y_poly]');
% T.Properties.VariableNames = {'In_x' 'B_y'};
% filename = 'latticeMatch.xlsx';
% writetable(T,filename)
