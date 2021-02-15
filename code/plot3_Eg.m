%% Read data from Excel
addpath(genpath('spectral_color_1'))
data = xlsread('../data/WZ_BInGaN_data_KG.xlsx','Master List - Full','C3:J115');

%% HSE Fit
B_y = data(1:40,1);
In_x = data(1:40,2);
Eg_HSE = data(1:40,8);

% Binary HSE Band Gaps [eV]
E_BN = 10.065;    % From Logan (for our functional) (direct G-G gap)
E_InN = 0.64;     % From Logan (for our functional)
E_GaN = 3.138;    % From Logan (for our functional)

% HSE fit parameters
tbl = table(In_x,B_y,Eg_HSE);
modelfunc = @(p,x) x(:,2).*E_BN + x(:,1).*E_InN + (1-x(:,1)-x(:,2)).*E_GaN + x(:,2).*x(:,1).*p(1) + x(:,1).*(1-x(:,1)-x(:,2)).*p(2) + x(:,2).*(1-x(:,1)-x(:,2)).*p(3);
p0 = [1 1 1];
model = fitnlm(tbl,modelfunc,p0);
p = model.Coefficients.Estimate
p_standardError = model.Coefficients.SE;
p_percentError = abs(p_standardError./p).*100

shift = 0.25;       %rigid shift to match experimental gap of GaN

%% Plot

% Create surface of best fit
figure(5)
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

%% Combined Figure
figure(72)
fs = 20;

X = X0;
Y = Y0;
Eg_fit0 = Y.*E_BN + X.*E_InN + (1-X-Y).*E_GaN + Y.*X.*p(1) + X.*(1-X-Y).*p(2) + Y.*(1-X-Y).*p(3);
subplot(2,1,2)
plot(X,Eg_fit0+shift,'LineWidth',2)
axis square
set(gca,'FontSize',fs,'Layer','top','XTick',[0:0.1:0.5],'YTick',[2.2:0.2:3.4],'YTickLabel',{'2.2','2.4','2.6','2.8','3.0','3.2',' '},'units','normalized','InnerPosition',[0 0.1100 1 0.4])
xlabel('Indium Mole Fraction','FontSize',fs)
ylabel('Band Gap [eV]','FontSize',fs)
xlim([0 0.5])
pos = get(gca,'position');

subplot('Position',[pos(1) pos(2)+pos(4) pos(3) pos(4)])
plot(X0,Y0,'LineWidth',2)
xlim([0 0.5])
ylim([0 0.5])
axis square
set(gca,'FontSize',fs,'XTick',[0:0.1:0.5],'XTickLabels',{},'YTickLabels',{' ','0.1','0.2','0.3','0.4','0.5'})
ylabel('Boron Mole Fraction','FontSize',fs)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
