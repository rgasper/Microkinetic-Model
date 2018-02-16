%% Batch reactor with a slightly more complicated mechanism
% (1) A + * <=> A*
% (2) A* + 2S* <=> B* + 2C*
% (3) B* <=> B + *
% (4) 2C* <=> C2 + 2*
% overall A <=> B + C2 catalyzed by S

clear all
close all

% alpha = (CtotS) / (Cmax * V), where CtotS (mol) total surface sites,
% Cmax (mol/m3) max fluid conc, V (m3) volume fluid
alpha = 0.1; % d'less, surface-to-gas capacity ratio
k1 = 0.1; % (m3/mol/s), forward rxn coeffic for step 1
km1 = 0.1; % (1/s), reverse rxn coeffic for step 1 (k-1)
k2 = 100.0; % (m6/mol2/s2)
km2 = 0.1; % (m6/mol2/s2)
k3 = 0.01; % (1/s)
km3 = 0.01; % (m3/mol/s)
k4 = 100.0; % (m3/mol/s)
km4 = 0.1; % (m6/mol2/s2)
Cmax = 1; % (mol/m3), max fluid-phase conc
k1p = k1*Cmax; % (1/s)
k2p = k2*(alpha*Cmax)^2; % (1/s) 
km2p = km2*(alpha*Cmax)^2; % (1/s)
km3p = km3*Cmax; % (1/s)
k4p = k4*(alpha*Cmax); % (1/s)
km4p = km4*Cmax*(alpha*Cmax); % (1/s)
% Keq = (k1*k2*k3*k4)/(km1*km2*km3*km4) is equilibrium constant for overall reaction
% Keq = [A]/2[C] isn't unitless as overall eqbm has unequal stoichiometry 
%, unitless version Keqp = Keq*Cmax
% equilibrium conversion should not change with Cmax since no change in total
% moles in overall reaction stoichiometry for this reaction
% set params array to send to user-written function in ode45 calls

params = [alpha k1p km1 k2p km2p k3 km3p k4p km4p];

%% integrate
%  y(1) = psi-A, y(2) = psi-B, y(3) =psi-C2
%  y(4) = theta-A, y(5) = theta-B, y(6) = theta-C
y0 = [1 0 0 0 0 0]; % initial conditions
tspan = [0 10000]; % time span
[t,y] = ode45('odes',tspan,y0,[],params);


%% Differentiate and Sensitivity Analysis
yp = diff(y,1);
tp = t(2:end,:);
yp = movmean(yp,25,1); %smooth it out

epsilon = 1e-6;
eq_cut = find_equilibrium(yp,epsilon);
%% Plot Reactor Progress
% no approximations applied, we have 3 ODEs
% and have no overall rate equation = function only of fluid-phase conc
%  dy1dt = alpha*( -k1p*y(1)*(1-y(3)-y(4)) + km1*y(3) );
%  dy2dt = alpha*( k3*y(4) - km3p*y(2)*(1-y(3)-y(4) );
%  dy3dt = k1p*y(1)*(1-y(3)-y(4)) - km1*y(3) - k2*y(3) + km2*y(4);
%  dy4dt = k2*y(3) - km2*y(4) - k3*y(4) + km3p*y(4)*(1-y(3)-y(4));
y1 = y(:,1);
y2 = y(:,2);
y3 = y(:,3);
y4 = y(:,4);
y5 = y(:,5);
y6 = y(:,6);
eqt = ones(length(t),1)*t(eq_cut);

figure(1)
plot(t,y1,'r', t,y2,'b', t,y3,'g',...
     t,y4,'r.',t,y5,'b.',t,y6,'g.')
tt = 'Reactor Detail';
title(tt)
xlabel('t')
ylim manual
hold on
ys = ylim;
eqy = linspace(ys(1),ys(2),length(eqt));
plot(eqt,eqy,'k--')
hold off
legend('\Psi_A','\Psi_B','\Psi_C',...
    '\theta_{A}','\theta_{B}','\theta_{C}','Eqbm.',... 
    'location','northEast')

%% Plot Derivatives
ypp = yp(101:end,:); %first hundred-ish throw off the plot detail
tpp = tp(101:end,:);
yp1 = ypp(:,1);
yp2 = ypp(:,2);
yp3 = ypp(:,3);
yp4 = ypp(:,4);
yp5 = ypp(:,5);
yp6 = ypp(:,6);
eqt = ones(length(tpp),1)*tpp(eq_cut-100);

figure(2)
plot(tpp,yp1,'r', tpp,yp2,'b', tpp,yp3,'g',...
     tpp,yp4,'r.',tpp,yp5,'b.',tpp,yp6,'g.')
tt = ...
 sprintf('Overall rates');
title(tt)
xlabel('t')
ylim manual
hold on
ys = ylim;
eqy = linspace(ys(1),ys(2),length(eqt));
plot(eqt,eqy,'k--')
hold off
legend(' \delta\Psi_A',' \delta\Psi_B',' \delta\Psi_C',...
    ' \delta\theta_{A}','\delta\theta_{B}','\delta\theta_{C}','Eqbm.',...
    'location','northEast')
