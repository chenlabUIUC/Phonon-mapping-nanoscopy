%% PMN on rhombic lattice
% PURPOSE:  Main code for mode analysis, plotting, and fitting of nearest
%           neighboring (NN) and angular (ANG) spring constants for rhombic
%           or square lattice.
%
% INPUT:    Drift corrected position matrix for rhombic or square lattice,
%           with shape of [N,2,t], N as number of nanoparticles, and t as
%           number of frames in 'Input_path'
%
% OUTPUT:   Figures and output data matrix in 'Plot_dir'
%           
% HISTORY:  Written by Chang Qian
% Last modified by Chang Qian on 07/09/2024

clear; close all; clc;
addpath('functions/')

%%
Input_path = 'Data drift rotation corrected/';
Plot_dir = 'Plots/';

activate_sets = [1:3];
% Input matrix name, tracking error (nm)
Input_params{1} = {'Cube_22mM.mat', 0.5};   % Estimated runtime 2-5 min
Input_params{2} = {'Cube_27mM.mat', 0.5};   % Estimated runtime 2-5 min
Input_params{3} = {'Cube_110mM.mat', 0.5};  % Estimated runtime <1 min

% Plotting parameters:
%   K space repetition
%   Max kx for 2D plot
%   Max frequency for 2D plot
%   Point size for 2D plot
Plotting_params{1} = {1, 0.1, 3.2, 50};
Plotting_params{2} = {1, 0.1, 1.7, 50};
Plotting_params{3} = {1, 0.1, 4.2, 50};

if ~exist(Plot_dir, 'dir'); mkdir(Plot_dir); end

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
    SegNT = 50;
    N_T2D = 200; % Require ~12 GB RAM
    [n1, n2] = deal(n2, n1);
    Q21 = pos;

    %% Lattice fitting, get Xref
    X = mean(Q21,3);
    E = @(x) LatticeError(X,x);
    x1 = fminsearch(E,x0);
    
    X1list=zeros(size(Q21,3),length(x1));
    for i=1:size(Q21,3) % lattice fitting frame by frame to find error in lattice constants.
        E2=@(x) LatticeError(Q21(:,:,i),x);
        X1list(i,:)=fminsearch(E2,x1);% save the lengths and angles of other lattice vector.
    end
    rot = atan(mean(X1list(:,2))/mean(X1list(:,1)));
    Q23=zeros(size(Q21));
    X1list3=X1list.*0.0;
    
    for i=1:size(Q21,3)
        Q23(:,1,i)=Q21(:,1,i)*cos(rot)+Q21(:,2,i)*sin(rot);
        Q23(:,2,i)=Q21(:,2,i)*cos(rot)-Q21(:,1,i)*sin(rot);
        X1list3(i,1)=X1list(i,1)*cos(rot)+X1list(i,2)*sin(rot);
        X1list3(i,2)=X1list(i,2)*cos(rot)-X1list(i,1)*sin(rot);
        X1list3(i,3)=X1list(i,3)*cos(rot)+X1list(i,4)*sin(rot);
        X1list3(i,4)=X1list(i,4)*cos(rot)-X1list(i,3)*sin(rot);
    end
    A1=mean(X1list3(:,1:2));
    A2=mean(X1list3(:,3:4));
    
    X = mean(Q23,3);
    Xref=zeros(length(X),2);
    N1=x1(5):x1(6);
    N2=x1(7):x1(8);
    for I=1:length(X)
        m=X(I,2)/A2(2);
        [~,b]=min((N2-m).^2);
        m=N2(b);
        n=(X(I,1)-m*A2(1))/A1(1);
        [~,b]=min((N1-n).^2);
        n=N1(b);
        Xref(I,1)=n*A1(1)+m*A2(1);
        Xref(I,2)=m*A2(2);
    end
    
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
    Path_rhomb = pathRhomb(b1, b2, SegN);
        % struct contain SegK_rec, SegDist_rec, SegPtDist, SegPtName
        % ready for interpolation
    
    % Box on plot
    Box_rhomb = boxRhomb(b1, b2);
    inBox = inpolygon(K(:,1), K(:,2), Box_rhomb(:,1), Box_rhomb(:,2));
    
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
        angle = acos(A1(1)*A2(1)/(A1mag*A2mag));
        l0 = (A1mag + A2mag)/2;
        
        % Rough estimation of spring constant
        x0 = [10^2,10^6];
        E = @(x) LatticeFitAng3(K,Wp,Wm,Wpperr_final',Wmmerr_final',Wpmerr,Echeck,x,angle,l0,Weight);
        
        Nspace1 = 100;
        Nspace2 = 100;
        x1 = logspace(-3,2,Nspace1);
        x2 = logspace(-2,8,Nspace2);
        [X1,X2] = meshgrid(x1,x2);
        X1 = transpose(X1);
        X2 = transpose(X2);
        Z = X1.*0;
        for i = 1:Nspace1
            if ~mod(i,10)
                disp(['Fitting: ',num2str(i),'/',num2str(Nspace1)])
            end
            for j = 1:Nspace2
                Z(i,j)=sqrt(E([X1(i,j),X2(i,j)]));
            end
        end
        
        Y2 = 2*X2./l0^2.*(2*angle-pi)^2;
        
        zcheck = islocalmin(Z).*islocalmin(Z,2);
        MinList = [X1(zcheck==1),X2(zcheck==1),Y2(zcheck==1),Z(zcheck==1)];
        X13 = X1(islocalmin(Z));
        X3 = X2(islocalmin(Z));
        Y3 = Y2(islocalmin(Z));
        Z2 = Z(islocalmin(Z));
        MinList = vertcat(MinList,[X13(1),X3(1),Y3(1),Z2(1)]);
        
        options = optimset('TolX',10^-12,'MaxIter',10000);
        S = size(MinList);
        
        MinList2 = zeros([S(1),4+3]);
        xm = x0;
        ErrorAng = sqrt(E(x0));
        for i = 1:S(1)
            x0 = [MinList(i,1),MinList(i,2)];
            
            % only consider Maxwellness<0.3 cases.
            if 2*x0(2)/l0^2*(2*angle-pi)^2/x0(1)>0.3; continue; end 
            
            x10 = fminsearch(E,x0,options);
            x10 = abs(x10);
            [d] = LatticeFitAng3check(K,Wp,Wm,Wpperr_final',Wmmerr_final',Wpmerr,Echeck,x10,angle,l0,Weight);
            ErrorAng2 = d; % sqrt(E(x10))
            MinList2(i,1) = x10(1); MinList2(i,2)=x10(2); MinList2(i,3)=2*x10(2)./l0^2.*(2*angle-pi)^2;
            MinList2(i,4) = d; %MinList2(i,5)=de; MinList2(i,6)=dm/d; MinList2(i,7)=dp/d;
            if ErrorAng2<ErrorAng
                ErrorAng = ErrorAng2;
                xm = x10;
            end
        end
        
        x10 = xm;
        x1 = x10;
        
        l0 = (A1mag+A2mag)/2;
        U2 = @(l,a) 2.*x1(1)./2.*(l-l0).^2+4.*x1(2).*(a-angle).^2.*(a+angle-pi).^2;
        kANGfit = [x1(1),x1(2)];
        [kANGS,kANGErr,kANGErrext] = LatticeEigAngErr3(K,Wp,Wm,Wpperr_final',Wmmerr_final',Wpmerr,Echeck,x1,angle,l0,0.001,Weight);
        kANGfitComp = [x1(1),2*x1(2)/l0^2*(2*angle-pi)^2];
        kANGfitCompErr = [kANGErr(1),2*kANGErr(2)/l0^2*(2*angle-pi)^2];
        kANGfitCompErrext = [kANGErrext(1),2*kANGErrext(2)/l0^2*(2*angle-pi)^2];
        kANGfitCompErrRSS = sqrt(kANGfitCompErr.^2+kANGfitCompErrext.^2);
        
        [lm,lp] = LatticeEigAng(x1,angle,0.5*(A1mag+A2mag));
        [Lpt,Lmt] = LatticePolAng(x1,angle,0.5*(A1mag+A2mag));
        
        % Modes, 2 bands.
        [KTx,KTy] = meshgrid(linspace(-1,1,N_T2D).*plot_klim, ...
                             linspace(-1,1,N_T2D).*plot_klim);
        KT = cat(2, reshape(KTx,[],1), reshape(KTy,[],1));
        Lp = sqrt(lp(KT(:,1),KT(:,2))); % lower
        Lm = sqrt(lm(KT(:,1),KT(:,2))); % upper
        Lp(Echeck==1) = 0.0;
        Lm(Echeck==1) = 0.0;
        
        % Polarization, 2 bands
        Lptlist=diag(Lmt(KT(:,1),KT(:,2)));
        Lmtlist=diag(Lpt(KT(:,1),KT(:,2)));
        Lptlist(Echeck==1)=0.0;
        Lmtlist(Echeck==1)=1.0;
        inBoxT = inpolygon(KT(:,1), KT(:,2), Box_rhomb(:,1), Box_rhomb(:,2));
        Lptlist(~inBoxT) = 0.5;
        Lmtlist(~inBoxT) = 0.5;
        
        %% 1D path interpolation and comparison
        unpackStruct(Path_rhomb);
        
        % Mode spectra, 1D curve, cubic interpolation
        % W_all = cat(2, Wp', dwp_mat, Wm', dwm_mat);
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
        PathT_rhomb = pathRhomb(b1, b2, SegNT);
        SegDistT = mean(PathT_rhomb.SegDist_rec,2);
        SegKT = PathT_rhomb.SegK_rec(:,:,1);
        
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
            plot(Box_rhomb(:,1),Box_rhomb(:,2),'k-')
            plot3(Box_rhomb(:,1),Box_rhomb(:,2),ones(size(Box_rhomb(:,1))).*plot_clim, 'k-', 'LineWidth',1)
        nexttile(2)
            standardModePlot(K, Wm, PLM, plot_klim, plot_clim, plot_pointSize, 'Experiment frequency, upper')
            plot(Box_rhomb(:,1),Box_rhomb(:,2),'k-')
            plot3(Box_rhomb(:,1),Box_rhomb(:,2),ones(size(Box_rhomb(:,1))).*plot_clim, 'k-', 'LineWidth',1)
        nexttile(3)
            % standardModePlot(KT, Lp, PLM, plot_klim, plot_clim, plot_pointSize/10, 'Model frequency, lower')
            imagesc(linspace(-1,1,N_T2D).*plot_klim, linspace(-1,1,N_T2D).*plot_klim, reshape(Lp,N_T2D,N_T2D));
            colormap(gca, PLM);
            title('Model frequency, lower')
            hold on;
            daspect([1,1,0.1]);
            axis([-1,1,-1,1].*plot_klim)
            set(gca,'YDir','normal')
            plot3(Box_rhomb(:,1),Box_rhomb(:,2),ones(size(Box_rhomb(:,1))).*plot_clim, 'k-', 'LineWidth',1)
        nexttile(4)
            % standardModePlot(KT, Lm, PLM, plot_klim, plot_clim, plot_pointSize/10, 'Model frequency, upper')
            imagesc(linspace(-1,1,N_T2D).*plot_klim, linspace(-1,1,N_T2D).*plot_klim, reshape(Lm,N_T2D,N_T2D));
            colormap(gca, PLM);
            title('Model frequency, upper')
            hold on;
            daspect([1,1,0.1]);
            axis([-1,1,-1,1].*plot_klim)
            set(gca,'YDir','normal')
            plot(Box_rhomb(:,1),Box_rhomb(:,2),'k-')
            plot3(Box_rhomb(:,1),Box_rhomb(:,2),ones(size(Box_rhomb(:,1))).*plot_clim, 'k-', 'LineWidth',1)
            colorbar;
        nexttile(5)
            standardModePlot(K(inBox,:), Lplist(inBox), BWR, plot_klim, 1, plot_pointSize, 'Experiment polarization, lower')
            plot(Box_rhomb(:,1),Box_rhomb(:,2),'k-')
            plot3(Box_rhomb(:,1),Box_rhomb(:,2),ones(size(Box_rhomb(:,1))).*plot_clim, 'k-', 'LineWidth',1)
        nexttile(6)
            standardModePlot(K(inBox,:), Lmlist(inBox), BWR, plot_klim, 1, plot_pointSize, 'Experiment polarization, upper')
            plot(Box_rhomb(:,1),Box_rhomb(:,2),'k-')
            plot3(Box_rhomb(:,1),Box_rhomb(:,2),ones(size(Box_rhomb(:,1))).*plot_clim, 'k-', 'LineWidth',1)
        nexttile(7)
            % standardModePlot(KT(inBoxT,:), Lptlist(inBoxT), BWR, plot_klim, 1, plot_pointSize/10, 'Model polarization, lower')
            imagesc(linspace(-1,1,N_T2D).*plot_klim, linspace(-1,1,N_T2D).*plot_klim, reshape(Lptlist,N_T2D,N_T2D));
            colormap(gca, BWR);
            caxis([0,1]);
            title('Model polarization, lower')
            hold on;
            daspect([1,1,0.1]);
            axis([-1,1,-1,1].*plot_klim)
            set(gca,'YDir','normal')
            plot(Box_rhomb(:,1),Box_rhomb(:,2),'k-')
            plot3(Box_rhomb(:,1),Box_rhomb(:,2),ones(size(Box_rhomb(:,1))).*plot_clim, 'k-', 'LineWidth',1)
        nexttile(8)
            % standardModePlot(KT(inBoxT,:), Lmtlist(inBoxT), BWR, plot_klim, 1, plot_pointSize/10, 'Model polarization, upper')
            imagesc(linspace(-1,1,N_T2D).*plot_klim, linspace(-1,1,N_T2D).*plot_klim, reshape(Lmtlist,N_T2D,N_T2D));
            colormap(gca, BWR);
            caxis([0,1]);
            title('Model polarization, upper')
            hold on;
            daspect([1,1,0.1]);
            axis([-1,1,-1,1].*plot_klim)
            set(gca,'YDir','normal')
            plot(Box_rhomb(:,1),Box_rhomb(:,2),'k-')
            plot3(Box_rhomb(:,1),Box_rhomb(:,2),ones(size(Box_rhomb(:,1))).*plot_clim, 'k-', 'LineWidth',1)
            colorbar;
        nexttile(9)
            cla;
            hold on; axis off;
            message = sprintf(['kNN: ',num2str(kANGfitComp(1)),' ± ',num2str(kANGfitCompErr(1)),'\n',...
                      'kANG: ',num2str(kANGfitComp(2)),' ± ',num2str(kANGfitCompErr(2))]);
            text(0,0,message)
        nexttile(10,[1,2])
            % hold on; box on;
            % errorbar(SegDistT,W_interp(:,1),W_interp(:,2),'vertical','o')
            % errorbar(SegDistT,W_interp(:,3),W_interp(:,4),'vertical','o')
            % plot(SegDistT,SegEigANG_lower,'-','LineWidth',2,'DisplayName','ANG model, lower')
            % plot(SegDistT,SegEigANG_upper,'-','LineWidth',2,'DisplayName','ANG model, upper')
            % ylabel('Frequency');
            % currentaxis = axis();currentaxis(3)=0;
            % axis(currentaxis);
            % xticks(SegPtDist); xticklabels(SegPtName); xlim([0,max(SegPtDist)]);
            % for i = 1:size(SegPtDist,1); xline(SegPtDist(i)); end
            % legend
            % legend({'Expt. mode lower','Expt. mode upper','ANG model lower','ANG model upper'},...
            %     'Location','southeast');

            hold on; box on;
            %  Frequency with errorbar
            xx = [mean(SegDist_rec,2); flip(mean(SegDist_rec,2))];
            
            errorbar(mean(SegDist_rec,2),W_interp(:,1),W_interp(:,2),'vertical','ko')
            
            yy = [W_interp(:,3)+W_interp(:,4) ; flip(W_interp(:,3)-W_interp(:,4)) ];
            fill(xx,yy,[1,1,1].*.6,'FaceAlpha',0.2,'LineStyle','none');
            plot(mean(SegDist_rec,2),W_interp(:,3), 'o-','Color',[1,1,1].*.4);
            %%%%%%%%%%%%%
            plot(SegDistT,SegEigANG_lower,'k-','LineWidth',2,'DisplayName','ANG model, lower')
            plot(SegDistT,SegEigANG_upper,'-','LineWidth',2,'DisplayName','ANG model, upper','Color',[1,1,1].*.4)
            ylabel('Frequency'); %ylim([-0.5e-3,4]); %axis auto
            currentaxis = axis();currentaxis(3)=0;
            axis(currentaxis);
            xticks(SegPtDist); xticklabels(SegPtName); xlim([0,max(SegPtDist)]);
            for i = 1:size(SegPtDist,1); xline(SegPtDist(i)); end
        
        saveFigures([Plot_dir, save_name, '_weight',num2str(Weight),'_full'])
        
        
        % Save all variables
        varNames = {'K','Wp','Wm','Lp','Lm','Lplist','Lmlist','Lptlist','Lmtlist','Wpperr','Wmmerr','Wpmerr','dwp_mat','dwm_mat',...
                    'SegDistT','W_interp','SegEigANG_lower','SegEigANG_upper','SegPtName','SegPtDist',...
                    'kANGfitComp','kANGfitCompErr','plot_klim','plot_clim','plot_pointSize','BWR','PLM',...
                    'Q23','U','Xref','A1','A2','b1','b2','n1','n2','Path_rhomb','Box_rhomb','Weight'};
        outputStruct = packStruct(varNames);
        save([Plot_dir, save_name, '_weight',num2str(Weight),'_mode.mat'],'outputStruct')
    
    end
end


%% Functions








