%% PMN on rhombic lattice

clear; close all; clc;
addpath('functions/');

%%
Input_path = 'Data drift rotation corrected/';
Plot_dir = 'Plots/';

activate_sets = [1];
% Input matrix name, tracking error (nm)
Input_params{1} = {'Rod.mat', 1};
Input_params{2} = {'Prism.mat', 5.7};

% Plotting parameters:
%   K space repetition
%   Max kx for 2D plot
%   Max frequency for 2D plot
%   Point size for 2D plot
Plotting_params{1} = {1, 0.12, 1.5, 50};
Plotting_params{2} = {1, 0.07, 0.4, 80};

if ~exist(Plot_dir, 'dir'); mkdir(Plot_dir); end

rot_mat = @(x) [cos(x), sin(x); -sin(x), cos(x)];

%%

for activate_set = activate_sets
    close all;

    %% Loading data
    In_name = Input_params{activate_set}{1};
    tracking_error = Input_params{activate_set}{2};
    load([Input_path, In_name])
    save_name = strrep(In_name, '.mat', '');

    plot_Krep = Plotting_params{activate_set}{1};
    plot_klim = Plotting_params{activate_set}{2};
    plot_clim = Plotting_params{activate_set}{3};
    plot_pointSize = Plotting_params{activate_set}{4};

    SegN = floor((n1+n2)/2);
    [n1, n2] = deal(n2, n1);
    Q21 = pos;

    %% Lattice fitting, get Xref
    [l_row, l_col] = getLabelsRowCol(Q21);
    col_offset = getOffset(n1);
    l_col = l_col - col_offset(l_row)';
    l_col = l_col - 1;
    l_row = l_row - 1;
    
    Qmean = mean(Q21,3);
    init_ind = find(ismember([l_row,l_col],[0,0],'rows'));
    a1_ind = find(ismember([l_row,l_col],[0,1],'rows'));
    a2_ind = find(ismember([l_row,l_col],[1,-1],'rows'));
    pt_init = Qmean(init_ind,:);
    pt_a1 = Qmean(a1_ind,:) - Qmean(init_ind,:);
    pt_a2 = Qmean(a2_ind,:) - Qmean(init_ind,:);
    lattice0 = [pt_init, pt_a1, pt_a2];
    
    E = @(x) LatticeError_hexa_new(x, Qmean, l_row, l_col);
    lattice_fit = fminsearch(E,lattice0);
    
    init = lattice_fit(1:2);
    a1 = lattice_fit(3:4);
    a2= lattice_fit(5:6);
    
    % Rotate everything, according to a1.
    rot = -AngleCalc([0,0],[1,0],a1)/180*pi;
    a1 = a1 * rot_mat(rot);
    a2 = a2 * rot_mat(rot);
    init = init * rot_mat(rot);
    Q23 = zeros(size(Q21));
    for i = 1:size(Q23,3)
        Q23(:,:,i) = Q21(:,:,i) * rot_mat(rot);
    end
    
    Q_lattice = l_row * a2 + l_col * a1 + init;
    
    A1 = a1;
    A2 = a2;
    Xref = Q_lattice;

    %% Preparation
    % U: displacement
    U = Q23 - mean(Q23,3); % Displacement from temporal averaged position.
    
    % K: resolvable reciprocal points
    [b1, b2] = reciprocalLattice(A1, A2);
        % Reciprocal lattice vector, b2 on y axis.
    
    K=zeros((2*plot_Krep*n1+1)*(2*plot_Krep*n2+1), 2);
    Echeck=zeros((2*plot_Krep*n1+1)*(2*plot_Krep*n2+1), 1);
        % Echeck is the Gamma points.
    ind=0;
    for i = -plot_Krep*n1 : plot_Krep*n1
        for j = -plot_Krep*n2 : plot_Krep*n2
            ind = ind + 1;
            K(ind,:) = i/n1*b1 + j/n2*b2;
            if mod(i,n1)==0 && mod(j,n2)==0
                Echeck(ind)=1;
            end
        end
    end
    
    % High symmetry path (consider symmetry)
    Path_hexa = pathHexa(b1, b2, SegN);
        % struct contain SegK_rec, SegDist_rec, SegPtDist, SegPtName
        % ready for interpolation
    
    % Box on plot
    Box_hexa = boxHexa(b1, b2);
    inBox = inpolygon(K(:,1),K(:,2),Box_hexa(:,1),Box_hexa(:,2));
    
    % Colormaps
    BWR = redblue(256);
    PLM = plasma(256);
    %% Equations
    Eqs = modeEquations(U, Xref);
    
    unpackStruct(Eqs);
    
    %% Calculate mode
    % Modes, 2 bands.
    Wp = wp(K(:,1),K(:,2)); % lower
    Wm = wm(K(:,1),K(:,2)); % upper
    Wp(Echeck==1) = 0.0;
    Wm(Echeck==1) = 0.0;
    
    % Error
    Wpperr = wpperr(K(:,1),K(:,2));
    Wmmerr = wmmerr(K(:,1),K(:,2));
    Wpmerr = wpmerr(K(:,1),K(:,2));
    Wpperr(Echeck==1) = 0.0;
    Wpmerr(Echeck==1) = 0.0;
    Wmmerr(Echeck==1) = 0.0;
    
    % Polarization, 2 bands
    Lplist = Lp(K(:,1),K(:,2));
    Lmlist = Lm(K(:,1),K(:,2));
    Lplist(Echeck==1) = 0.0;
    Lmlist(Echeck==1) = 1.0;
    
    %% New frequency error
    [dwp_mat, dwm_mat, other] = errorPosition_20240606(U, Xref, K, Wp, Wm, M, tracking_error);
    
    unpackStruct(other);
    
    Wpperr_final = sqrt(Wpperr'.^2 + dwp_mat.^2);
    Wmmerr_final = sqrt(Wmmerr'.^2 + dwm_mat.^2);
    
    %% Model fitting
    for Weight = [1] % Weighted fitting to both phonon branches
        A1mag = sqrt(sum(A1.^2));
        A2mag = sqrt(sum(A2.^2));
        angle = acos(A1(1)*A2(1) / (A1mag*A2mag));
        l0 = (A1mag + A2mag) / 2;
        E=@(x) LatticeFit_hexa(K,Wp,Wm,Wpperr_final',Wmmerr_final',Wpmerr,Echeck,x,l0,Weight);
        
        % Rough fit: spring constant scan
        Nspace1=500;
        x1=logspace(-5,5,Nspace1);
        Z=x1.*0;
        for i=1:Nspace1
            Z(i)=E(x1(i));
        end
        Z = real(Z);
        zcheck = islocalmin(Z);
        xlist = x1(zcheck)';
        
        % Fine fit: fminsearch
        options = optimset('TolX',10^-12,'MaxIter',10000);
        for i=1:size(xlist,1)
            x10=fminsearch(E,xlist(i,1),options);
            xlist(i,2:3) = [x10,E(x10)];
        end
        
        Emin = min(xlist(:,3));
        kfit = xlist(xlist(:,3) == Emin,2);
        U2 = @(l,a) 2.*x1(1)./2.*(l-l0).^2+4.*x1(2).*(a-angle).^2.*(a+angle-pi).^2;
        % [kANGS,kfitErr,kANGErrext]=LatticeEigAngErr_hexa(K,Wp,Wm,Wpperr,Wmmerr,Wpmerr,Echeck,kfit,l0,0.001,Weight);
        [kANGS,kfitErr,kANGErrext]=LatticeEigAngErr_hexa_20240630(K,Wp,Wm,Wpperr_final',Wmmerr_final',Wpmerr,Echeck,kfit,l0,0.001,0);
        
        [lp,lm]=LatticeEig_hexa(kfit,l0);
        [Lpt,Lmt]=LatticePol_hexa(kfit,l0);
        
        % Modes, 2 bands.
        Lp=sqrt(lp(K(:,1),K(:,2))); % lower
        Lm=sqrt(lm(K(:,1),K(:,2))); % upper
        Lp(Echeck==1)=0.0;
        Lm(Echeck==1)=0.0;
        
        % Polarization, 2 bands
        Lptlist=diag(Lmt(K(:,1),K(:,2)));
        Lmtlist=diag(Lpt(K(:,1),K(:,2)));
        Lptlist(Echeck==1)=0.0;
        Lmtlist(Echeck==1)=1.0;
        
        %% 1D path interpolation and comparison
        unpackStruct(Path_hexa);
        
        % Mode spectra, 1D curve, cubic interpolation
        W_all = cat(2, Wp', Wpperr_final, Wm', Wmmerr_final);
        W_interp = [];
        for i = 1:size(W_all,2)
            temp = [];
            for j = 1:size(SegK_rec,3)
                temp = cat(2,temp,griddata(K(:,1), K(:, 2), W_all(:,i), SegK_rec(:,1,j), SegK_rec(:,2,j),'cubic'));
            end
            W_interp = cat(2,W_interp,mean(temp,2));
        end
        
        % ANG model
        SegDistT = mean(SegDist_rec,2);
        SegKT = SegK_rec(:,:,1);
        
        SegEigANG_lower = zeros(size(SegDistT));
        SegEigANG_upper = zeros(size(SegDistT));
        for i=1:size(SegKT,1)
            SegEigANG_lower(i)=sqrt(lm(SegKT(i,1),SegKT(i,2)));
            SegEigANG_upper(i)=sqrt(lp(SegKT(i,1),SegKT(i,2)));
        end

        %% Output
        % Big mode plot, 4*3
        figure(); tiledlayout(3,4,"TileSpacing",'compact','Padding','compact');
        set(gcf,'Position', [30 30 1600 900])
        nexttile(1)
            standardModePlot(K, Wp, PLM, plot_klim, plot_clim, plot_pointSize, 'Experiment frequency, lower')
            plot(Box_hexa(:,1),Box_hexa(:,2),'k-')
        nexttile(2)
            standardModePlot(K, Wm, PLM, plot_klim, plot_clim, plot_pointSize, 'Experiment frequency, upper')
            plot(Box_hexa(:,1),Box_hexa(:,2),'k-')
        nexttile(3)
            standardModePlot(K, Lm, PLM, plot_klim, plot_clim, plot_pointSize, 'Model frequency, lower')
            plot(Box_hexa(:,1),Box_hexa(:,2),'k-')
        nexttile(4)
            standardModePlot(K, Lp, PLM, plot_klim, plot_clim, plot_pointSize, 'Model frequency, upper')
            plot(Box_hexa(:,1),Box_hexa(:,2),'k-')
            colorbar;
        nexttile(5)
            standardModePlot(K(inBox,:), Lplist(inBox), BWR, plot_klim, 1, plot_pointSize, 'Experiment polarization, lower')
            plot(Box_hexa(:,1),Box_hexa(:,2),'k-')
        nexttile(6)
            standardModePlot(K(inBox,:), Lmlist(inBox), BWR, plot_klim, 1, plot_pointSize, 'Experiment polarization, upper')
            plot(Box_hexa(:,1),Box_hexa(:,2),'k-')
        nexttile(7)
            standardModePlot(K(inBox,:), Lptlist(inBox).*0, BWR, plot_klim, 1, plot_pointSize, 'Model polarization, lower')
            plot(Box_hexa(:,1),Box_hexa(:,2),'k-')
        nexttile(8)
            standardModePlot(K(inBox,:), Lmtlist(inBox).*0+1, BWR, plot_klim, 1, plot_pointSize, 'Model polarization, upper')
            plot(Box_hexa(:,1),Box_hexa(:,2),'k-')
            colorbar;
        nexttile(9)
            cla;
            hold on; axis off;
            message = sprintf(['kFit: ',num2str(kfit), ' ± ', num2str(kfitErr)]);
            text(0,0,message)
        nexttile(10,[1,2])
            hold on; box on
            %%%%%%%%%%%%%
            %  Frequency with errorbar
            xx = [SegDistT; flip(SegDistT)];
            
            errorbar(SegDistT,W_interp(:,1),W_interp(:,2),'vertical','ko')
            
            yy = [W_interp(:,3)+W_interp(:,4) ; flip(W_interp(:,3)-W_interp(:,4)) ];
            fill(xx,yy,[1,1,1].*.6,'FaceAlpha',0.2,'LineStyle','none');
            plot(SegDistT,W_interp(:,3), 'o-','Color',[1,1,1].*.4);
            %%%%%%%%%%%%%
            plot(SegDistT,SegEigANG_lower,'k-','LineWidth',2,'DisplayName','ANG model, lower')
            plot(SegDistT,SegEigANG_upper,'-','LineWidth',2,'DisplayName','ANG model, upper','Color',[1,1,1].*.4)
            ylabel('Frequency'); %ylim([-0.5e-3,4]); %axis auto
            currentaxis = axis();currentaxis(3)=0;
            axis(currentaxis);
            xticks(SegPtDist); xticklabels(SegPtName); xlim([0,max(SegPtDist)]);
            for i = 1:size(SegPtDist,1); xline(SegPtDist(i)); end
            legend
            legend({'Expt. mode lower','','Expt. mode upper','ANG model lower','ANG model upper'},...
                'Location','southeast');
        
        saveFigures([Plot_dir, save_name, '_weight',num2str(Weight), '_full'])
        
        
        % Save all variables
        varNames = {'K','Wp','Wm','Lp','Lm','Lplist','Lmlist','Lptlist','Lmtlist','Wpperr','Wmmerr','Wpmerr','dwp_mat','dwm_mat',...
                    'SegDistT','W_interp','SegEigANG_lower','SegEigANG_upper','SegPtName','SegPtDist',...
                    'kfit','kfitErr','plot_klim','plot_clim','plot_pointSize','BWR','PLM',...
                    'Q23','U','Xref','A1','A2','b1','b2','n1','n2','Path_hexa','Box_hexa','Weight'};
        outputStruct = packStruct(varNames);
        save([Plot_dir, save_name, '_weight',num2str(Weight), '_mode.mat'],'outputStruct')

    end
end





%% Functions








