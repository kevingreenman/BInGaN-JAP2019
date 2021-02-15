data = xlsread('../data/WZ_BInGaN_data_KG.xlsx','Master List - Full','C3:F115');
B_y = data(:,1);
In_x = data(:,2);
c = data(:,4);

% Lattice Parameters [A]
c_BN = 4.22136;
c_InN = 5.75288;
c_GaN = 5.207065026;

% Get fit parameters
tbl = table(In_x,B_y,c);
modelfunc = @(p,x) x(:,2).*c_BN + x(:,1).*c_InN + (1-x(:,1)-x(:,2)).*c_GaN + x(:,2).*x(:,1).*p(1) + x(:,1).*(1-x(:,1)-x(:,2)).*p(2) + x(:,2).*(1-x(:,1)-x(:,2)).*p(3); %bowing model
p0 = [1 1 1];
model = fitnlm(tbl,modelfunc,p0);
p = model.Coefficients.Estimate
p_standardError = model.Coefficients.SE;
p_percentError = abs(p_standardError./p).*100

% Surfaces
x = linspace(0,0.7,500);
y = x;
[X,Y] = meshgrid(x,y);
c_fit = Y.*c_BN + X.*c_InN + (1-X-Y).*c_GaN + Y.*X.*p(1) + X.*(1-X-Y).*p(2) + Y.*(1-X-Y).*p(3);
c_fit(Y>1-X) = NaN;     %Show only physical region
cGaN = c_GaN*ones(length(X));
cGaN(Y>1-X) = NaN;     %Show only physical region

relErr = (c_fit-cGaN)./cGaN;
figure(32)
surf(X,Y,relErr);
xlabel('Indium Mole Fraction')
ylabel('Boron Mole Fraction')
zlabel('(c - c_{GaN})/c_{GaN}')
title('Lattice Parameter (c): Relative Difference with GaN')
axis equal
xlim([0 0.7])
ylim([0 0.5])
shading flat
grid off
box on
c = colorbar;
set(c,'FontSize',15)
hold on
view(0,90) %show plot in 2D
set(gca,'FontSize',15,'Layer','top')
colorcet('D1');
caxis([-0.04 0.04])
set(c,'YTick',[-0.04:0.01:0.04],'TickLabels',{'-4%', '-3%', '-2%', '-1%', ' 0%', ' 1%', ' 2%', ' 3%', ' 4%'},'TickLength',0)


%Isobars
spacing = 0.01;
for iso = -0.15:spacing:0.05 %For all isobars, use -0.1:spacing:0.1
    X_iso = X(abs(relErr - iso) < 1E-4);
    Y_iso = Y(abs(relErr - iso) < 1E-4);
    if ~isempty(X_iso) && ~isempty(Y_iso)
        X_iso_interp = linspace(0,0.7,500);
        Y_iso_interp = interp1(X_iso,Y_iso,X_iso_interp);
        plot3(X_iso_interp,Y_iso_interp,ones(size(X_iso_interp)).*max(max(relErr))+0.001,'--k')
    end
end

%Plot black line along where relative error is zero
relErr0 = abs(relErr) < 5E-5;
X0 = X(relErr0);
Y0 = Y(relErr0);
plot3(X0,Y0,zeros(size(X0)),'k')

set(findobj(gca, 'Type', 'Line', 'Linestyle', '-'), 'LineWidth', 2); %Make line thicker

%White Cover-Up
[X,Y] = meshgrid(linspace(0,0.7,1000),linspace(0,0.7,1000));
Z = ones(size(X)).*max(max(relErr))+0.001;
Z(Y<X & Y>0.5*X & Y<1-X) = NaN;
q = surf(X,Y,Z,'FaceColor','w');
set(q,'edgecolor','none')
