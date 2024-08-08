%% Fit ANG Spring
% PURPOSE:  This code shows how we fit angular spring constant from coarse
%           grained (CG) energy calculation.
%
% INPUT:    Results from CG energy calculation.
%           Energy_Minima.mat: 3-column matrix 'M' of [I, bl, ag] as ionic
%           strength [mM], equilibrium bond length [nm] and angle [deg].
%           ANGSpring.mat: 3-column matrix 'rec' of [I, ag, E] as ionic
%           strength [mM], angle [deg], and energy [kbT]. Bond length is
%           fixed at l_m, which is the equilibrium bond length at the
%           corresponding ionic strength condition.
%
% OUTPUT:   fitResult_ANG matrix:
%           5-column matrix of [I, k_ANG, dk_ANG, ag_min, ag_max] as ionic
%           strength [mM], angular spring constant [kbT/nm^2], fitting
%           error (95% confidence), fitting range (min and max) [deg].
%           
% HISTORY:  Written by Chang Qian
% Last modified by Chang Qian on 09/19/2023
%%
clear; close all; clc;

I0 = [22,27,110]; % mM
load('Data/Cube_Energy_Minima.mat')
load('Data/Cube_ANGSpring.mat')

figure(1); hold on;
fitResults_ANG = [];

for I = I0
    select = rec(:,1)==single(I);
    t = rec(select,[2,3]);
    Em = min(t(:,2));
    select = t(:,2)<Em+1;
    t2 = t(select,:);
    BL0 = M(find(I0==I),2);
    AG0 = M(find(I0==I),3);
    disp([min(t2(:,1)),max(t2(:,1))]);
    
    xfit = [72:0.02:108];
    [kANG,x0,yfit,rsq,se] = fitANGSpringConst((t2(:,1)-AG0).^2.*(t2(:,1)+AG0-180).^2,...
                                t2(:,2),AG0,xfit);
    
    kANG = kANG  /4*2/BL0^2*(2*AG0/180*pi-pi)^2;
    se(1) = se(1)/4*2/BL0^2*(2*AG0/180*pi-pi)^2;
    
    figure(1);
        plot(t(:,1), t(:,2),'o-');
        plot(t2(:,1),t2(:,2),'k-','Linewidth',2);
        plot(xfit, yfit, 'k--');
        
    fitResults_ANG = [fitResults_ANG; I, kANG, se(1),...
                              min(t2(:,1)),max(t2(:,1))]; 
                          
end
ylim([-4,4])


function [kANG,x0,yfit,rsq,se] = fitANGSpringConst(L,E,AG0,xrange_fit)

    select = 1:length(L);
    x_temp = L(select)*(1/180*pi)^4;
    E_temp = E(select);
    x_temp = double(x_temp); E_temp = double(E_temp);
    w = exp(-E_temp);
    
    [xData, yData, weights] = prepareCurveData( x_temp, E_temp, w );
        ft = fittype( 'poly1' );
        opts = fitoptions( 'Method', 'LinearLeastSquares' );
        opts.Weights = weights;
    [fitresult, gof] = fit( xData, yData, ft, opts );
    
    a = fitresult.p1;
    b = fitresult.p2;
    Rsquare = gof.rsquare;
    
    alpha = 0.95; df = length(x_temp);
    ci = confint(fitresult, alpha);
    t = tinv((1+alpha)/2, df); 
	se = (ci(2,:)-ci(1,:)) ./ (2*t); % Standard Error
    
    g = @(x) a.*(x-AG0).^2.*(x+AG0-180).^2*(1/180*pi).^4 + b;
    
    kANG = a;
    x0 = b*180/pi; if x0>90; x0 = 180-x0; end
    yfit = g(xrange_fit);
    rsq = Rsquare;
    
end









