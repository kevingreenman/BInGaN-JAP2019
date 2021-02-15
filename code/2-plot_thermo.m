%% Setup
data = xlsread('../data/WZ_BInGaN_data_KG.xlsx','Master List - Full','C3:I87'); %Wurtzite
B_y = data(:,1);
In_x = data(:,2);
H = data(:,6);
T = data(:,7);

x = linspace(0,0.7,500);
y = x;
[X,Y] = meshgrid(x,y);

% Calculate entropy and then critical temp at every point
Z = 1-X-Y;
S = ones(size(X));
for i=1:size(X,1)
    for j=1:size(X,2)
        if Y(i,j)>1-X(i,j)
            S(i,j) = NaN;
        else
            if X(i,j)==0
                S_x = 0;
            else
                S_x = X(i,j).*log(X(i,j));
            end
            if Y(i,j)==0
                S_y = 0;
            else
                S_y = Y(i,j).*log(Y(i,j));
            end
            if Z(i,j)==0
                S_z = 0;
            else
                S_z = Z(i,j).*log(Z(i,j));
            end
            
            S(i,j) = S_x + S_y + S_z;
        end
    end
end


%% Enthalpy of Mixing Per Cation

% Get fit parameters
tbl = table(In_x,B_y,H);
modelfunc = @(p,x) p(1).*x(:,1).*x(:,2) + p(2).*x(:,1).*(1-x(:,1)-x(:,2)) + p(3).*x(:,2).*(1-x(:,1)-x(:,2));  % regular solution model
p0 = [1 1 1];
model = fitnlm(tbl,modelfunc,p0);
p = model.Coefficients.Estimate
p_standardError = model.Coefficients.SE;
p_percentError = abs(p_standardError./p).*100

% Plot surface
figure(11)
H_fit = p(1).*X.*Y + p(2).*X.*(1-X-Y) + p(3).*Y.*(1-X-Y);
H_fit(Y>1-X) = NaN;     %Show only physical region
surf(X,Y,H_fit,'EdgeColor','none');
hold on
[C,h] = contour3(X,Y,H_fit,0.05:0.1:0.65,'--k');
contour3(X,Y,H_fit+1,1.05:0.1:1.65,'--k');
xlabel('Indium Mole Fraction')
ylabel('Boron Mole Fraction')
title('Enthalpy of Mixing [eV/cation]')
hold on
colorcet('L9')
caxis([0 0.7])
colorbar
view([0 90])
set(gca, 'Layer', 'top');
box on
axis equal
xlim([0 0.7])
ylim([0 0.5])
grid off

% %White Cover-Up       %need to comment out this section if running entire script at once
[X,Y] = meshgrid(linspace(0,0.7,1000),linspace(0,0.7,1000));
Z = ones(size(X)).*max(max(H_fit))+0.001;
Z(Y<X & Y>0.5*X & Y<1-X) = NaN;
q = surf(X,Y,Z,'FaceColor','w');
set(q,'edgecolor','none')

set(gca,'FontSize',15,'Layer','top')

%% Phase-Transition Temperature

x = linspace(0,0.7,500);
y = x;
[X,Y] = meshgrid(x,y);
T_fit = (H_fit.*300)./(-0.025.*S);

% Plot surface
figure(12)
surf(X,Y,T_fit,'EdgeColor','none');
hold on
contour3(X,Y,T_fit,[2000:1000:12000],'--k');
contour3(X,Y,T_fit+12000,[14000:1000:24000],'--k');
xlabel('Indium Mole Fraction')
ylabel('Boron Mole Fraction')
title('Phase-Transition Temperature [K]')
colorcet('L9')
colorbar
view([0 90])
set(gca, 'Layer', 'top');
box on
axis equal
xlim([0 0.7])
ylim([0 0.5])
caxis([1000 12000])
grid off

%White Cover-Up
[X,Y] = meshgrid(linspace(0,0.7,1000),linspace(0,0.7,1000));
Z = ones(size(X)).*max(max(T_fit))+0.001;
Z(Y<X & Y>0.5*X & Y<1-X) = NaN;
q = surf(X,Y,Z,'FaceColor','w');
set(q,'edgecolor','none')

set(gca,'FontSize',15,'Layer','top')

% Calculate entropy at every point
Z0 = 1-X0-Y0;
S0 = ones(length(X0),1);
for i=1:length(X0)
    if Y0(i)>1-X0(i)
        S0(i) = NaN;
    else
        if X0(i)==0
            S_x = 0;
        else
            S_x = X0(i).*log(X0(i));
        end
        if Y0(i)==0
            S_y = 0;
        else
            S_y = Y0(i).*log(Y0(i));
        end
        if Z0(i)==0
            S_z = 0;
        else
            S_z = Z0(i).*log(Z0(i));
        end
        
        S0(i) = S_x + S_y + S_z;
    end
end

%% combined figure
fs = 20;
figure(71)
H_fit0 = p(1).*X0.*Y0 + p(2).*Y0.*(1-X0-Y0) + p(3).*X0.*(1-X0-Y0);
subplot(2,1,2)
plot(X0,H_fit0,'LineWidth',2)
axis square
set(gca,'FontSize',fs,'Layer','top','XTick',[0:0.1:0.5],'YTick',[0:0.1:0.7],'units','normalized','InnerPosition',[0 0.1100 1 0.4],'YTickLabels',{'0','0.1','0.2','0.3','0.4','0.5','  0.6',' '})
xlabel('Indium Mole Fraction','FontSize',fs)
ylabel('Enthalpy of Mixing [eV/cation]','FontSize',fs)
xlim([0 0.5])
ylim([0 0.7])
pos = get(gca,'position');

T_fit0 = (H_fit0.*300)./(-0.025.*S0);
subplot('Position',[pos(1) pos(2)+pos(4) pos(3) pos(4)])
plot(X0(3:end),T_fit0(3:end),'LineWidth',2)
axis square
set(gca,'FontSize',fs,'Layer','top','XTick',[0:0.1:0.5],'YTick',[1000:2000:11000],'XTickLabels',{},'YTickLabels',{' ','3.0','5.0','7.0','9.0','11.0'})
xlim([0 0.5])
ylim([1000 11000])
ylabel('Temperature [10^3 K]','FontSize',fs)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);


