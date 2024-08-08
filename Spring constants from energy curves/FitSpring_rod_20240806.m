%% Spring constant fitting, Rod
% PURPOSE:  This code shows how we fit spring constant from coarse grained
%           (CG) energy calculation.
%
% INPUT:    Results from CG energy calculation.
%           Data/Rod.mat: include variables 'dList' and 'E_sums', as
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
load('Data/Rod.mat')
x_curve = dList;
y_curve = E_sums';

%%
% Select range:
%   Left: define by 1 kBT energy
%   Right: define by displacement
    Eu = 1; % kBT 
    r_displacement = 1.2; % nm

Em = min(y_curve);
rmin = x_curve(find(y_curve==min(y_curve),1));
select = y_curve<=Em+Eu;
xtemp = x_curve(select);
xmin = min(xtemp); % Left bound
xmax = rmin + r_displacement; % Right bound

select = x_curve >= xmin & x_curve <= xmax; % Fitting range of the curve

% Curve fitting
xtemp = x_curve(select);
ytemp = y_curve(select);
w = exp(-ytemp);
[fitresult, gof] = createFit(xtemp, ytemp, w, rmin, Em);
kfit = fitresult.a;

% Visualizationn
figure()
set(gcf,'Position',[150 150 300 300])
plot(x_curve, y_curve);
hold on;
plot(x_curve(select), y_curve(select),'k-','LineWidth',2)
a = axis();
xx = linspace(a(1), a(2),200);
plot(xx, fitresult(xx),'r--')
axis(a);
title('Interaction curve')
xlabel('r (nm)')
ylabel('U (k_BT)')


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

