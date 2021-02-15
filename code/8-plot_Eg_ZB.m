%% Read data from Excel
addpath(genpath('spectral_color_1'))
data = xlsread('../data/ZB_BInGaN_data_KG.xlsx','Master List - Full','C3:I89');

%% HSE Fit
B_y = data(1:35,1);
In_x = data(1:35,2);
Eg_HSE = data(1:35,7);

% Binary HSE Band Gaps [eV]
% Experimental
% E_BN = 6.25;      % ZB: http://www.ioffe.ru/SVA/NSM/Semicond/BN/basic.html#zinc%20blende
% E_InN = 0.7;      % ZB: https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=7488132
% E_GaN = 3.24;     % ZB: http://www.semiconductors.co.uk/nitrides.htm#GaN
% Calculated
E_BN = 10.549028;    % From Logan (for our functional)
E_InN = 0.446521;   % From Logan (for our functional)
E_GaN = 2.927436;   % From Logan (for our functional)

% HSE fit parameters
tbl = table(In_x,B_y,Eg_HSE);
modelfunc = @(p,x) x(:,2).*E_BN + x(:,1).*E_InN + (1-x(:,1)-x(:,2)).*E_GaN + x(:,2).*x(:,1).*p(1) + x(:,1).*(1-x(:,1)-x(:,2)).*p(2) + x(:,2).*(1-x(:,1)-x(:,2)).*p(3);
p0 = [1 1 1];
model = fitnlm(tbl,modelfunc,p0);
p = model.Coefficients.Estimate
p_standardError = model.Coefficients.SE;
p_percentError = abs(p_standardError./p).*100

shift = 0.31;       %rigid shift to match experimental gap of GaN

%% Plot

% Create surface of best fit
figure(6)
set(gcf,'Position',[100 100 700 450])
[X,Y] = meshgrid(linspace(0.001,0.699,1000),linspace(0.001,0.499,1000));
Eg_fit = Y.*E_BN + X.*E_InN + (1-X-Y).*E_GaN + Y.*X.*p(1) + X.*(1-X-Y).*p(2) + Y.*(1-X-Y).*p(3);
X = reshape(X,[],1);
Y = reshape(Y,[],1);
Eg_fit = reshape(Eg_fit,[],1)+shift;
colors = spectrumRGB(1240./Eg_fit);
colors = permute(colors,[2 3 1]);

% Plot surface of best fit in 2D
scatter(X,Y,Eg_fit,colors,'filled')
xlabel('Indium Mole Fraction')
ylabel('Boron Mole Fraction')
title('Calculated Band Gap [eV]')
hold on
view(0,90) %show plot in 2D
set(gca,'FontSize',15,'Layer','top')
set(gca,'OuterPosition',[0 0 0.8 1])
set(gca, 'Layer', 'top');
box on
axis equal

%Isobars
spacing = 0.2;
for iso = 2.0:spacing:3.0
    X_iso = X(abs(Eg_fit - iso) < 5E-4);
    Y_iso = Y(abs(Eg_fit - iso) < 5E-4);
    [Ysorted, Yorder] = sort(Y_iso);
    Xsorted = X_iso(Yorder,:);
    if iso == 2.8 || iso==3.0
        plot3(Xsorted,Ysorted,zeros(size(X_iso))+0.1,'--k')     % was white in colored region
    else
        plot3(Xsorted,Ysorted,zeros(size(X_iso))+0.1,'--k')
    end
end

%White Cover-Up
[X,Y] = meshgrid(linspace(0,0.7,1000),linspace(0,0.7,1000));
Z = ones(size(X)).*max(max(relErr))+0.001;
Z(Y<X & Y>0.5*X) = NaN;
q = surf(X,Y,Z,'FaceColor','w');
set(q,'edgecolor','none')
Z2 = ones(size(X)).*max(max(relErr))+0.101;
Z2(X+Y<1) = NaN;
q2 = surf(X,Y,Z2,'FaceColor','w');
set(q2,'edgecolor','none')

xlim([0 0.7])
ylim([0 0.5])
spectrumLabel(gca) %insert colorbar
