%% Fit NN Spring
% PURPOSE:  This code shows how we fit nearest neighbour (NN) spring   
%           constants from coarse grained (CG) energy calculation.
%
% INPUT:    Results from CG energy calculation.
%           Energy_Minima.mat: 3-column matrix 'M' of [I, bl, ag] as ionic
%           strength [mM], equilibrium bond length [nm] and angle [deg].
%           NNSpring.mat: 3-column matrix 'rec' of [I, bl, E] as ionic
%           strength [mM], bond length [nm], and energy [kbT]. Angle is
%           fixed at theta_m, which is the equilibrium angle at the
%           corresponding ionic strength condition.
%
% OUTPUT:   fitResult_ANG matrix:
%           5-column matrix of [I, k_NN, dk_NN, bl_min, bl_max] as ionic
%           strength [nm], angular spring constant [kbT/nm^2], fitting
%           error (95% confidence), fitting range (min and max) [nm].
%           
% HISTORY:  Written by Chang Qian
% Last modified by Chang Qian on 09/19/2023
%%
clear; close all; clc;

I0 = [22,27,110]; % mM
load('Data/Cube_Energy_Minima.mat')
load('Data/Cube_NNSpring.mat')

figure(1); hold on;
fitResults_NN = [];

for I = I0
    select = rec(:,1)==single(I);
    t = rec(select,[2,3]);
    Em = min(t(:,2));
    select = t(:,2)<Em+1;
    t2 = t(select,:);
    bl0 = M(find(I0==I),2);
    disp([min(t2(:,1)),max(t2(:,1))]);
    
    xfit = [68:0.02:82];
    [kNN,x0,yfit,rsq,se] = fitNNSpringConst((t2(:,1)-bl0).^2,t2(:,2),bl0,xfit);
    
    figure(1);
        plot(t(:,1), t(:,2),'o-');
        plot(t2(:,1),t2(:,2),'k-','Linewidth',2);
        plot(xfit, yfit, 'k--');
        
    fitResults_NN = [fitResults_NN; I, kNN*2, se(1)*2,...
                              min(t2(:,1)),max(t2(:,1))];
                          
end
ylim([-22,13])


function [kNN,x0,yfit,rsq,se] = fitNNSpringConst(L,E,bl0,xrange_fit)

    select = 1:length(L);
    x_temp = L(select);
    E_temp = E(select);
    x_temp = double(x_temp); E_temp = double(E_temp);
    w = exp(-E_temp);
    
    [xData, yData, weights] = prepareCurveData(x_temp, E_temp, w );
        ft = fittype( 'poly1' );
        opts = fitoptions( 'Method', 'LinearLeastSquares' );
        opts.Weights = weights;
    [fitresult, gof] = fit( xData, yData, ft, opts );
    
    p1 = fitresult.p1;
    p2 = fitresult.p2;
    Rsquare = gof.rsquare;
    
    alpha = 0.95; df = length(x_temp);
    ci = confint(fitresult, alpha);
    t = tinv((1+alpha)/2, df); 
	se = (ci(2,:)-ci(1,:)) ./ (2*t); % Standard Error
    
    g = @(x) p1.*(x-bl0).^2 + p2;
    
    kNN = p1;
    x0 = p2/2/p1;
    yfit = g(xrange_fit);
    rsq = Rsquare;
    
end









