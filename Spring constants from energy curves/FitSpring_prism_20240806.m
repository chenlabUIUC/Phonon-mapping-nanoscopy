%% Spring constant fitting, Prism
% PURPOSE:  This code shows how we fit spring constant from coarse grained
%           (CG) energy calculation.
%
% INPUT:    Results from CG energy calculation.
%           Data/Prism.mat: include variables 'dList' and 'E_sums', as
%           vectors of interparticle center-to-center distance and
%           interaction energy as a sum of van der Waals and electrostatic
%           interactions.
%
% OUTPUT:   Fitting result: spring constant in variable 'kfit', and plot of
%           original energy curve (blue), fitting range (black), and fitted
%           energy curve (red).
%           
% HISTORY:  Written by Chang Qian
% Last modified by Chang Qian on 08/08/2024

clear; close all; clc;

%% Load data
load('Data/Prism.mat')
x_curve = dList;
y_curve = E_sums;

%%
% Select range:
%   Right: Define by displacement
%   Left: Same energy level as right bound
    r_displacement = 5.8;

rmin = x_curve(find(y_curve==min(y_curve),1));
xmax = rmin + r_displacement; % Right bound
select = x_curve <= xmax;
xtemp = x_curve(select);
ytemp = y_curve(select);

Em = min(y_curve);
EM = ytemp(end);
select = y_curve <= EM; % Left bound
xtemp = x_curve(select);
ytemp = y_curve(select);
w = exp(-ytemp);

[fitresult, gof] = createFit(xtemp', ytemp, w, rmin, Em);
kfit = fitresult.a;

figure()
set(gcf,'Position',[150 150 500 500])
nexttile()
plot(x_curve, y_curve);
hold on;
plot(x_curve(select), y_curve(select),'k-','LineWidth',2)
a = axis();
xx = linspace(a(1), a(2),200);
plot(xx, fitresult(xx),'r-')
axis(a);
title('Interaction curve')
xlabel('r (nm)')
ylabel('U (kBT)')
xlim([110,140])
ylim([-0.35,0.4])


%% Functions

function [fitresult, gof] = createFit(xtemp, ytemp, w, xmin, Emin)
    [xData, yData, weights] = prepareCurveData( xtemp, ytemp, w );

    % Set up fittype and options.
    ft = fittype( ['a/2*(x-',num2str(xmin),')^2+c'], 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.StartPoint = [0.933993247757551 Emin];
    opts.Weights = weights;
    
    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft, opts );
end

