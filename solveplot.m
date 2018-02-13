%% Batch reactor with simple mechanism
% (1) A + * <=> A*
% (2) A* <=> B*
% (3) B* <=> C + *
% overall A <=> C catalyzed by S
clear all
close all
% alpha = (CtotS) / (Cmax * V), where CtotS (mol) total surface sites,
% Cmax (mol/m3) max fluid conc, V (m3) volume fluid
alpha = 0.01; % d'less, surface-to-gas capacity ratio
k1 = 1.0; % (m3/mol/s), forward rxn coeffic for step 1
km1 = 0.1; % (1/s), reverse rxn coeffic for step 1 (k-1)
k2 = 1.0; % (1/s)
km2 = 0.1; % (m3/mol/s)
k3 = 0.1; % (1/s)
km3 = 0.1; % (m3/mol/s)
Cmax = 1; % (mol/m3), max fluid-phase conc
k1p = k1*Cmax; % (1/s)
km3p = km3*Cmax; % (1/s)
% Keq = (k1*k2*k3)/(km1*km2*km3) is equilibrium constant for overall reaction
% equilibrium conversion should not change with Cmax since no change in total
% moles in overall reaction stoichiometry for this reaction
% set params array to send to user-written function in ode45 calls
params = [alpha k1p km1 k2 km2 k3 km3p];
% integrate
y0 = [1 0 0 0]; % initial conditions
tspan = [0 4000]; % time span
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
eqt = ones(length(t),1)*t(eq_cut);

figure(1)
plot(t,y1,'b',t,y2,'r',t,y3,'k--',t,y4,'m--')
tt = 'Reactor Detail';
title(tt)
xlabel('t')
ylim manual
hold on
ys = ylim;
eqy = linspace(ys(1),ys(2),length(eqt));
plot(eqt,eqy,'g--')
hold off
legend('\Psi_A','\Psi_C',' \theta_{A}','\theta_{B}','Eqbm.',... 
    'location','northEast')

%% Plot Derivatives
ypp = yp(101:end,:); %first hundred-ish throw off the plot detail
tpp = tp(101:end,:);
yp1 = ypp(:,1);
yp2 = ypp(:,2);
yp3 = ypp(:,3);
yp4 = ypp(:,4);
eqt = ones(length(tpp),1)*tpp(eq_cut);

figure(2)
plot(tpp,yp1,'b',tpp,yp2,'r',tpp,yp3,'k',tpp,yp4,'m')
tt = ...
 sprintf('Overall rates');
title(tt)
xlabel('t')
ylim manual
hold on
ys = ylim;
eqy = linspace(ys(1),ys(2),length(eqt));
plot(eqt,eqy,'g--')
hold off
legend(' \delta\Psi_A',' \delta\Psi_C',' \delta\theta_{A}',...
    '\delta\theta_{B}','Eqbm.','location','northEast')
