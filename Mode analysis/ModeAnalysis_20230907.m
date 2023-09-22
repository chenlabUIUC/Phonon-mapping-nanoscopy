%% Mode analysis.m
% PURPOSE:  Perform mode analysis to tracked lattice region. Mode
%           calculation work for all lattices in principle, but lattice
%           fitting and mode fitting work for rhombic lattice only.
%
% INCLUDE:  Evolutionaryoptomize.m: evolutionary optimization to find the
%           best lattice parameters.
%           framecompare3.m: loss for drift/rotation correction.
%           LatticeEigAng.m: calculate eigenvalue of discrete mechanical
%           model.
%           LatticeEigAngErr3.m: calculate error of spring constants in
%           discrete mechanical model.
%           LatticeError.m: loss function for lattice fitting
%           LatticeFitAng3.m, LatticeFitAng3check.m: loss function for
%           discrete mechanical moodel fitting.
%           LatticePolAng.m: calculate polarization of discrete mechanical
%           model.
%           plasma.m: generate plasma colormap.
%           redblue.m: generate red-white-blue colormap.
%
% INPUT:    Position matrix of tracked particles, with shape (n,2,t), with n
%           as the number of particles, 2 as x and y direction, t as total number of
%           frames.
%
% OUTPUT:   Figure similar to Fig.2 in the manuscript.
%
% HISTORY:  Written by Chang Qian， Ethan Stanifer
% Last modified by Chang Qian on 09/19/2023
%%
clear; close all; clc;
addpath('functions/')

%%%%%%%%%%%%%%%%%% Input %%%%%%%%%%%%%%%%%%%%%%%%%%%
Input_name = 'Pos_22mM.mat';

n1 = 7; % columns
n2 = 13; % rows

P = 1; % Repeatition of mode plot in k space

% Plot drift correction/lattice fitting results
check_intermediate_step = 1;
% Plot mode
plot_mode = 1;

% Lattice fitting result, not changing results but just accelerating
% calculation time for this case. Comment this for the first run.
x0 = [78.1319948285825 6.81198309532632e-05 17.2078956589178 74.3959852609863 ...
      -2.99999948693110 3.00000051306890 -5.99999948496878 6.00000051503122];
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load data
if ~exist('Plots/','dir'); mkdir('Plots/'); end
load(Input_name);

Q = Pos;
L = size(Q,3);
Num = size(Q,1);

%% Drift correction
Q20 = Q;
for i = 1:L; Q20(:,:,i) = Q(:,:,i) - mean(Q(:,:,i)); end

%% Rotation correction
R0 = zeros(1,L); % Initial guess of rotation for each frame
C = @(x) framecompare3(x,Q20);
options = optimset('Display','iter','PlotFcns',@optimplotfval);
R = fminsearch(C,R0,options);
R(1) = 0; % Rotation of the 1st frame force to be zero

% Get position after rotation.
Q21=zeros(size(Q20));
for i=1:L
    Q21(:,1,i)=Q20(:,1,i)*cos(R(i))+Q20(:,2,i)*sin(R(i));
    Q21(:,2,i)=Q20(:,2,i)*cos(R(i))-Q20(:,1,i)*sin(R(i));
end

if plot_mode
    f = figure();
    f.Position(3) = f.Position(4)*2;
    subplot(1,2,1)
    plot(reshape(Q(:,1,:),Num,[])',reshape(Q(:,2,:),Num,[])','k')
    axis equal
    title('Trajectories')
    subplot(1,2,2)
    plot(reshape(Q21(:,1,:),Num,[])',reshape(Q21(:,2,:),Num,[])','k')
    axis equal
    title('Drift corrected')
    drawnow;
end


%% Lattice fitting
%mean position
X = mean(Q21,3);
E = @(x) LatticeError(X,x);

if ~exist('x0','var')
    s1=(n1-1)/2;
    s2=(n2-1)/2;
    x0=[100,0,0,100,-s1,s1,-s2,s2];
    lb=[-200,-200,-200,-200,-6,3,-6,3];
    ub=[200,200,200,200,-3,6,-3,6];
    mutation=[10,10,10,10,0.1,0.1,0.1,0.1];
    % Largest time investment on exploring the parameter space: length and
    % direction of lattice vectors, number of repeating units along the two
    % vectors.
    x = Evolutionaryoptomize(E,lb,ub,x0,2000,300,1600,100,mutation,4000);
else
    x = x0; 
end
x1 = fminsearch(E,x);

% Rotate lattice to align lattice vector a1 on +x axis.
X1list=zeros(L,length(x1));
for i=1:L % lattice fitting frame by frame to find error in lattice constants.
    E2=@(x) LatticeError(Q21(:,:,i),x);
    X1list(i,:)=fminsearch(E2,x1);% save the lengths and angles of other lattice vector.
end
rot = atan(mean(X1list(:,2))/mean(X1list(:,1)));
Q23=zeros(size(Q21));
X1list3=X1list.*0.0;
%rotate such that a1 is in the x direction on average
%two options: rot - A1y=0, angle1 - average a1 angle is 0
%we choose rot.
for i=1:L
    Q23(:,1,i)=Q21(:,1,i)*cos(rot)+Q21(:,2,i)*sin(rot);
    Q23(:,2,i)=Q21(:,2,i)*cos(rot)-Q21(:,1,i)*sin(rot);
    X1list3(i,1)=X1list(i,1)*cos(rot)+X1list(i,2)*sin(rot);
    X1list3(i,2)=X1list(i,2)*cos(rot)-X1list(i,1)*sin(rot);
    X1list3(i,3)=X1list(i,3)*cos(rot)+X1list(i,4)*sin(rot);
    X1list3(i,4)=X1list(i,4)*cos(rot)-X1list(i,3)*sin(rot);
end

% Parameters of the lattice
A1=mean(X1list3(:,1:2));
A2=mean(X1list3(:,3:4));
A1s=std(X1list3(:,1:2));
A2s=std(X1list3(:,3:4));
A1mag=sqrt(A1(1)^2+A1(2)^2);
A2mag=sqrt(A2(1)^2+A2(2)^2);

% Parameters for the reciprocal lattice
b1=2*pi*[A2(2),-A2(1)]./A2(2)./A1mag;
b2=2*pi*[0,1]./A2(2);

% Echeck: Gamma points in k space.
K=zeros([(2*P*n1+1)*(2*P*n2+1),2]);
Echeck=zeros([(2*P*n1+1)*(2*P*n2+1),1]);
I=0;
for n=-P*n1:P*n1
    for m=-P*n2:P*n2
        I=I+1;
        K(I,1)=n/n1*b1(1)+m/n2*b2(1);
        K(I,2)=n/n1*b1(2)+m/n2*b2(2);
        if and(mod(n,n1)==0,mod(m,n2)==0)
            Echeck(I)=1;
        end
    end
end

% Preparation for high symmetry path
BoxX=0.5.*[b1(1)+b2(1),b1(1)-b2(1),-b1(1)-b2(1),-b1(1)+b2(1),b1(1)+b2(1)];
BoxY=0.5.*[b1(2)+b2(2),b1(2)-b2(2),-b1(2)-b2(2),-b1(2)+b2(2),b1(2)+b2(2)];
inRhomb = inpolygon(K(:,1),K(:,2),BoxX,BoxY);

% References lattice (ideal lattice)
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

% Variation (displacement from averaged position)
U=Q23-mean(Q23,3);

%% Mode Calculation
% Equations
Imag=0+1j;

P = size(U);
uq=@(k1,k2,a) 1/sqrt(Num)*transpose(exp(-1.0*Imag*(k1.*transpose(Xref(:,1))+k2.*transpose(Xref(:,2))))*reshape(U(:,a,:),[P(1),P(3)]));

M=@(k1,k2,a,b) sum(uq(-k1,-k2,a).*uq(k1,k2,b))/L;
At = @(k1,k2) uq(-k1,-k2,1).*uq(k1,k2,1);
Bt = @(k1,k2) real(uq(-k1,-k2,1).*uq(k1,k2,2));
Ct = @(k1,k2) imag(uq(-k1,-k2,1).*uq(k1,k2,2));
Dt = @(k1,k2) uq(-k1,-k2,2).*uq(k1,k2,2);

SQ= @(k1,k2) sqrt(((M(k1,k2,1,1)+M(k1,k2,2,2))./2).^2-(M(k1,k2,1,1).*M(k1,k2,2,2)-M(k1,k2,1,2).*M(k1,k2,2,1)));

% Frequency of two bands
wp= @(k1,k2) (((M(k1,k2,1,1)+M(k1,k2,2,2))./2)+SQ(k1,k2)).^(-1/2);
wm= @(k1,k2) (((M(k1,k2,1,1)+M(k1,k2,2,2))./2)-SQ(k1,k2)).^(-1/2);

% Error of frequency
CvAA = @(k1,k2) mean(At(k1,k2).*At(k1,k2))-mean(At(k1,k2)).*mean(At(k1,k2));
CvAB = @(k1,k2) mean(At(k1,k2).*Bt(k1,k2))-mean(At(k1,k2)).*mean(Bt(k1,k2));
CvAC = @(k1,k2) mean(At(k1,k2).*Ct(k1,k2))-mean(At(k1,k2)).*mean(Ct(k1,k2));
CvAD = @(k1,k2) mean(At(k1,k2).*Dt(k1,k2))-mean(At(k1,k2)).*mean(Dt(k1,k2));
CvBB = @(k1,k2) mean(Bt(k1,k2).*Bt(k1,k2))-mean(Bt(k1,k2)).*mean(Bt(k1,k2));
CvBC = @(k1,k2) mean(Bt(k1,k2).*Ct(k1,k2))-mean(Bt(k1,k2)).*mean(Ct(k1,k2));
CvBD = @(k1,k2) mean(Bt(k1,k2).*Dt(k1,k2))-mean(Bt(k1,k2)).*mean(Dt(k1,k2));
CvCC = @(k1,k2) mean(Ct(k1,k2).*Ct(k1,k2))-mean(Ct(k1,k2)).*mean(Ct(k1,k2));
CvCD = @(k1,k2) mean(Ct(k1,k2).*Dt(k1,k2))-mean(Ct(k1,k2)).*mean(Dt(k1,k2));
CvDD = @(k1,k2) mean(Dt(k1,k2).*Dt(k1,k2))-mean(Dt(k1,k2)).*mean(Dt(k1,k2));

dpA=@(k1,k2) (-0.25*(wp(k1,k2)).^3).*(1+0.5*(M(k1,k2,1,1)-M(k1,k2,2,2))./SQ(k1,k2));
dmA=@(k1,k2) (-0.25*(wm(k1,k2)).^3).*(1-0.5*(M(k1,k2,1,1)-M(k1,k2,2,2))./SQ(k1,k2));
dpD=@(k1,k2) (-0.25*(wp(k1,k2)).^3).*(1-0.5*(M(k1,k2,1,1)-M(k1,k2,2,2))./SQ(k1,k2));
dmD=@(k1,k2) (-0.25*(wm(k1,k2)).^3).*(1+0.5*(M(k1,k2,1,1)-M(k1,k2,2,2))./SQ(k1,k2));

dpB=@(k1,k2) (-0.25*(wp(k1,k2)).^3).*(2*real(M(k1,k2,1,2))./SQ(k1,k2));
dmB=@(k1,k2) (-0.25*(wm(k1,k2)).^3).*(2*real(M(k1,k2,1,2))./SQ(k1,k2));
dpC=@(k1,k2) (-0.25*(wp(k1,k2)).^3).*(2*imag(M(k1,k2,1,2))./SQ(k1,k2));
dmC=@(k1,k2) (-0.25*(wm(k1,k2)).^3).*(2*imag(M(k1,k2,1,2))./SQ(k1,k2));

wpperr=@(k1,k2) wp(k1,k2).^2 ./ 2 ./ (L-1);
wmmerr=@(k1,k2) wm(k1,k2).^2 ./ 2 ./ (L-1);
wpmerr=@(k1,k2) wm(k1,k2) .* 0;

% Polarization:

Ppx= @(k1,k2) M(k1,k2,1,2)./sqrt(norm(M(k1,k2,1,1)-wp(k1,k2).^2).^2+norm(M(k1,k2,1,2)).^2);
Ppy= @(k1,k2) -1.0*(M(k1,k2,1,1)-wp(k1,k2).^2)./sqrt(norm(M(k1,k2,1,1)-wp(k1,k2).^2).^2+norm(M(k1,k2,1,2)).^2);
Pmx= @(k1,k2) M(k1,k2,1,2)./sqrt(norm(M(k1,k2,1,1)-wm(k1,k2).^2).^2+norm(M(k1,k2,1,2)).^2);
Pmy= @(k1,k2) -1.0*(M(k1,k2,1,1)-wm(k1,k2).^2)./sqrt(norm(M(k1,k2,1,1)-wm(k1,k2).^2).^2+norm(M(k1,k2,1,2)).^2);

% Compute:
Wp=wp(K(:,1),K(:,2));
Wm=wm(K(:,1),K(:,2));
Wpperr=wpperr(K(:,1),K(:,2));
Wmmerr=wmmerr(K(:,1),K(:,2));
Wpmerr=wpmerr(K(:,1),K(:,2));

PpxList=Ppx(K(:,1),K(:,2));
PpyList=Ppy(K(:,1),K(:,2));
PmxList=Pmx(K(:,1),K(:,2));
PmyList=Pmy(K(:,1),K(:,2));

Wp(Echeck==1)=0.0;
Wm(Echeck==1)=0.0;
Wp2=Wp;
Wm2=Wm;

yp= @(k1,k2) (((M(k1,k2,2,2)-M(k1,k2,1,1))./2)+SQ(k1,k2))./sqrt(M(k1,k2,1,2).*M(k1,k2,2,1)+(((M(k1,k2,2,2)-M(k1,k2,1,1))./2)+SQ(k1,k2)).^2);
ym= @(k1,k2) (((M(k1,k2,2,2)-M(k1,k2,1,1))./2)-SQ(k1,k2))./sqrt(M(k1,k2,1,2).*M(k1,k2,2,1)+(((M(k1,k2,2,2)-M(k1,k2,1,1))./2)-SQ(k1,k2)).^2);
xp= @(k1,k2) M(k1,k2,1,2)./sqrt(M(k1,k2,1,2).*M(k1,k2,2,1)+(((M(k1,k2,2,2)-M(k1,k2,1,1))./2)+SQ(k1,k2)).^2);
xm= @(k1,k2) M(k1,k2,1,2)./sqrt(M(k1,k2,1,2).*M(k1,k2,2,1)+(((M(k1,k2,2,2)-M(k1,k2,1,1))./2)-SQ(k1,k2)).^2);

Lm= @(k1,k2) (xm(k1,k2).*transpose(k1)+ym(k1,k2).*transpose(k2))...
          .*conj((xm(k1,k2).*transpose(k1)+ym(k1,k2).*transpose(k2)))...
          ./transpose(k1.^2+k2.^2);
Lp= @(k1,k2) (xp(k1,k2).*transpose(k1)+yp(k1,k2).*transpose(k2))...
          .*conj((xp(k1,k2).*transpose(k1)+yp(k1,k2).*transpose(k2)))...
          ./transpose(k1.^2+k2.^2);

Lplist=Lp(K(:,1),K(:,2));
Lmlist=Lm(K(:,1),K(:,2));
Lplist(Echeck==1)=0.0;
Lmlist(Echeck==1)=1.0;

BWR = redblue(256);
PLM = plasma(256);

if plot_mode
    f = figure(3);
    f.Position(3) = f.Position(4)*2;
    % Compare corrected position and reference lattice
    subplot(1,2,1); hold on;
    for i=1:L; plot(Q23(:,1,i),Q23(:,2,i),'.g'); end
    daspect([1 1 1]); box on;axis equal; drawnow;
    plot(Xref(:,1),Xref(:,2),'.r')
    title('Corrected position vs. ideal lattice')
    % Distribution of variation
    subplot(1,2,2); hold on;
    for i=1:L;plot(U(:,1,i),U(:,2,i),'.b');end;axis equal; drawnow;
    title('Variation')
    
    f = figure(4);
    f.Position = [70 70 1600 900];
    % Frequency, lower
    subplot(3,4,1);
%     set(gca,'Position',[0.025 0.55, 0.225,0.4]);
    scatter3(K(:,1),K(:,2),Wp2,50,Wp2,'filled'); 
    view(2)
    daspect([1 1 1]); box on; axis equal; drawnow;
    xlim([-0.1,0.1]); ylim([-0.1,0.1]);
    colormap(gca,PLM)
    colorbar;
    caxis([0,2.2])
    title('Frequency, lower')
    % Frequency, upper
    subplot(3,4,2);
%     set(gca,'Position',[0.275 0.55, 0.225,0.4]);
    scatter3(K(:,1),K(:,2),Wm2,50,Wm2,'filled');
    view(2)
    daspect([1 1 1]); box on;axis equal; drawnow;
    xlim([-0.1,0.1]); ylim([-0.1,0.1]);
    colormap(gca,PLM)
    colorbar;
    caxis([0,2.2])
    title('Frequency, upper')
    % Polarization, lower
    subplot(3,4,5); 
%     set(gca,'Position',[0.025 0.05, 0.225,0.4]);
%     scatter3(K(:,1),K(:,2),Lplist,50,Lplist,'filled'); 
    scatter3(K(inRhomb,1),K(inRhomb,2),Lplist(inRhomb),50,Lplist(inRhomb),'filled'); 
    view(2)
    hold on; plot(BoxX,BoxY,'k');
    daspect([1 1 1]); box on;axis equal; drawnow;
    xlim([-0.1,0.1]); ylim([-0.1,0.1]);
    colorbar
    caxis([0,1])
    colormap(gca,BWR)
    title('Polarization, lower')
    % Polarization, upper
    subplot(3,4,6); 
%     set(gca,'Position',[0.275 0.05, 0.225,0.4]);
%     scatter3(K(:,1),K(:,2),Lmlist,50,Lmlist,'filled'); 
    scatter3(K(inRhomb,1),K(inRhomb,2),Lmlist(inRhomb),50,Lmlist(inRhomb),'filled'); 
    view(2)
    hold on; plot(BoxX,BoxY,'k');
    daspect([1 1 1]); box on;axis equal; drawnow;
    xlim([-0.1,0.1]); ylim([-0.1,0.1]);
    colorbar
    caxis([0,1])
    colormap(gca,BWR)
    title('Polarization, upper')
    
    % Polarization, upper+lower, should be 1
    if check_intermediate_step
        figure();
        scatter3(K(inRhomb,1),K(inRhomb,2),Lmlist(inRhomb)+Lplist(inRhomb),...
                    50,Lmlist(inRhomb)+Lplist(inRhomb),'filled'); 
        view(2)
        hold on; plot(BoxX,BoxY,'k');
        daspect([1 1 1]); box on;axis equal; 
        colorbar
        colormap(BWR);
        title('Polarization, sum')
        drawnow;
    end
end

%% Discrete mechanical model fitting
Weight=0;
angle=acos(A1(1)*A2(1)/(A1mag*A2mag));
l0=(A1mag+A2mag)/2;

% Rough estimation of spring constant
x0=[10^2,10^6];
E=@(x) LatticeFitAng3(K,Wp,Wm,Wpperr,Wmmerr,Wpmerr,Echeck,x,angle,l0,Weight);

Nspace1=200;
Nspace2=200;
x1=logspace(-3,2,Nspace1);
x2=logspace(-2,8,Nspace2);
[X1,X2]=meshgrid(x1,x2);
X1=transpose(X1);
X2=transpose(X2);
Z=X1.*0;
for i=1:Nspace1
    disp(['Fitting: ',num2str(i),'/',num2str(Nspace1)])
    for j=1:Nspace2
        Z(i,j)=sqrt(E([X1(i,j),X2(i,j)]));
    end
end

Y2=2*X2./l0^2.*(2*angle-pi)^2;

zcheck=islocalmin(Z).*islocalmin(Z,2);
MinList=[X1(zcheck==1),X2(zcheck==1),Y2(zcheck==1),Z(zcheck==1)];
X13=X1(islocalmin(Z));
X3=X2(islocalmin(Z));
Y3=Y2(islocalmin(Z));
Z2=Z(islocalmin(Z));
MinList=vertcat(MinList,[X13(1),X3(1),Y3(1),Z2(1)]);

if check_intermediate_step
    figure(); hold on;
    surf(X1,Y2,Z,'FaceAlpha',0.75,'EdgeColor','none')
    set(gca, 'YScale', 'log')
    set(gca, 'XScale', 'log')
    set(gca, 'ZScale', 'log')
    set(gca,'ColorScale','log')
    view(3);
    title('Least Squares Residual')
    xlabel('kNN (kbT/nm^2)') 
    ylabel('kANG (kbT/nm^2)') 
    zlabel('Residual') 
    plot3(X1(islocalmin(Z)),Y2(islocalmin(Z)),Z(islocalmin(Z)),'.k');
    plot3(X1(islocalmin(Z,2)),Y2(islocalmin(Z,2)),Z(islocalmin(Z,2)),'.b');
    plot3(X1(zcheck==1),Y2(zcheck==1),Z(zcheck==1),'*r');
end

options=optimset('TolX',10^-12,'MaxIter',10000);
S = size(MinList);

MinList2=zeros([S(1),4+3]);

xm=x0;
ErrorAng=sqrt(E(x0));
for i=1:S(1)
    x0=[MinList(i,1),MinList(i,2)];
    
    % only consider Maxwellness<0.3 cases.
    if 2*x0(2)/l0^2*(2*angle-pi)^2/x0(1)>0.3; continue; end 
    
    x10=fminsearch(E,x0,options);
    x10=abs(x10);
    [d,dm,dp,de]=LatticeFitAng3check(K,Wp,Wm,Wpperr,Wmmerr,Wpmerr,Echeck,x10,angle,l0,Weight);
    ErrorAng2=d; % sqrt(E(x10))
    MinList2(i,1)=x10(1); MinList2(i,2)=x10(2); MinList2(i,3)=2*x10(2)./l0^2.*(2*angle-pi)^2;
    MinList2(i,4)=d; MinList2(i,5)=de; MinList2(i,6)=dm/d; MinList2(i,7)=dp/d;
    if ErrorAng2<ErrorAng
        ErrorAng=ErrorAng2;
        xm=x10;
    end
end

x10=xm;
x1=x10;

l0=(A1mag+A2mag)/2;
U2 = @(l,a) 2.*x1(1)./2.*(l-l0).^2+4.*x1(2).*(a-angle).^2.*(a+angle-pi).^2;
kANGfit=[x1(1),x1(2)];
[kANGS,kANGErr,kANGErrext]=LatticeEigAngErr3(K,Wp,Wm,Wpperr,Wmmerr,Wpmerr,Echeck,x1,angle,l0,0.001,Weight);
kANGfitComp=[x1(1),2*x1(2)/l0^2*(2*angle-pi)^2];
kANGfitCompErr=[kANGErr(1),2*kANGErr(2)/l0^2*(2*angle-pi)^2];
kANGfitCompErrext=[kANGErrext(1),2*kANGErrext(2)/l0^2*(2*angle-pi)^2];
kANGfitCompErrRSS=sqrt(kANGfitCompErr.^2+kANGfitCompErrext.^2);
Llist=linspace(l0-20,l0+20,2000);
Alist=linspace(pi/3,2*pi/3,1000);
Ulist=U2(Llist,transpose(Alist));
Ulist2=Ulist;
Ulist2(Ulist2>4*x1(2).*((pi/2-angle).^2.*(pi/2+angle-pi).^2)*2.5)=1/0.0;

EnergyBarrierAngular=4*x1(2).*((pi/2-angle).^2.*(pi/2+angle-pi).^2);
[lp,lm]=LatticeEigAng(x1,angle,0.5*(A1mag+A2mag));
[Lpt,Lmt]=LatticePolAng(x1,angle,0.5*(A1mag+A2mag));

Lptlist=diag(Lmt(K(:,1),K(:,2)));
Lmtlist=diag(Lpt(K(:,1),K(:,2)));
Lptlist(Echeck==1)=0.0;
Lmtlist(Echeck==1)=1.0;

if plot_mode
    figure(4);
    % Model frequency, upper
    subplot(3,4,4);
%     set(gca,'Position',[0.775 0.55, 0.225,0.4]);
    hold on; 
    E0=lp(K(:,1),K(:,2));
    scatter3(K(:,1),K(:,2),sqrt(E0),50,sqrt(E0),'filled'); 
    view(2)
    daspect([1 1 1]); box on;axis equal; drawnow;
    colormap(gca,PLM)
    colorbar
    caxis([0,2.2])
    xlim([-0.1,0.1]); ylim([-0.1,0.1]);
    title('Model Frequency, upper')
    % Model frequency, lower
    subplot(3,4,3);
%     set(gca,'Position',[0.525 0.55, 0.225,0.4]);
    hold on; 
    E0=lm(K(:,1),K(:,2));
    scatter3(K(:,1),K(:,2),sqrt(E0),50,sqrt(E0),'filled'); 
    view(2)
    daspect([1 1 1]); box on;axis equal; drawnow;
    colormap(gca,PLM)
    colorbar
    caxis([0,2.2])
    xlim([-0.1,0.1]); ylim([-0.1,0.1]);
    title('Model frequency, lower')
    % Model polarization, lower
    subplot(3,4,7);
%     set(gca,'Position',[0.525 0.05, 0.225,0.4]);
    scatter3(K(inRhomb,1),K(inRhomb,2),Lptlist(inRhomb),50,Lptlist(inRhomb),'filled'); 
    view(2)
    hold on; plot(BoxX,BoxY,'k');
    daspect([1 1 1]); box on;axis equal; drawnow;
    colorbar
    caxis([0,1])
    colormap(gca,BWR)
    xlim([-0.1,0.1]); ylim([-0.1,0.1]);
    title('Model polarization, lower')
    % Model polarization, upper
    subplot(3,4,8);
%     set(gca,'Position',[0.775 0.05, 0.225,0.4]);
    scatter3(K(inRhomb,1),K(inRhomb,2),Lmtlist(inRhomb),50,Lmtlist(inRhomb),'filled'); 
    view(2)
    hold on; plot(BoxX,BoxY,'k');
    daspect([1 1 1]); box on;axis equal; drawnow;
    colorbar
    caxis([0,1])
    colormap(gca,BWR)
    xlim([-0.1,0.1]); ylim([-0.1,0.1]);
    title('Model frequency, upper')
    if check_intermediate_step
        % Model polarization, upper+lower, should be 1
        figure();
        scatter3(K(inRhomb,1),K(inRhomb,2),Lmtlist(inRhomb)+Lptlist(inRhomb),...
                    50,Lmtlist(inRhomb)+Lptlist(inRhomb),'filled'); 
        view(2)
        hold on; plot(BoxX,BoxY,'k');
        daspect([1 1 1]); box on;axis equal; 
        colorbar
        %caxis([0,1])
        colormap(gca,BWR);drawnow;
        title('Model polarization, sum')
    end
end

%% Compare expt. vs. model
    % Path: Gamma-->R1-->M-->Gamma-->R2-->M
    SegN = 10;
    Gamma = [0,0]; M = b1/2; R1 = (b1+b2)/2; R2 = (b1-b2)/2;
    R1mag = sqrt(sum(R1.^2)); R2mag = sqrt(sum(R2.^2));
    if R1mag>R2mag; R3=R1; R1=R2; R2=R3; end
    R1mag = sqrt(sum(R1.^2)); R2mag = sqrt(sum(R2.^2));
    Mmag = sqrt(sum(M.^2)); R1Mmag = sqrt(sum((R1-M).^2)); R2Mmag = sqrt(sum((R2-M).^2));
    SegPt = [Gamma;R1;M;Gamma;R2;M]; % Main.
    for i = 1:size(SegPt,1)
        if i==1
            SegK = [0,0];
            SegDist = [0];
        else
            kx=linspace(real(SegPt(i-1,1)),real(SegPt(i,1)),SegN);
            ky=linspace(real(SegPt(i-1,2)),real(SegPt(i,2)),SegN);
            kxy = [kx;ky]';
            SegK = cat(1,SegK,kxy(2:end,:));
            Dist = sqrt(sum((kxy-SegPt(i-1,:)).^2,2))+SegDist(end);
            SegDist = cat(1,SegDist,Dist(2:end));
        end
    end
    SegPtDist = SegDist(1:SegN-1:end);
    SegPtName = {'\Gamma','R_1','M','\Gamma','R_2','M'};


    normVec = [R1(2),-R1(1)]; normVecMag = sqrt(sum(normVec.^2));
    normVec = normVec/normVecMag;
    ReflectMat = [1-2*normVec(1)^2,-2*normVec(1)*normVec(2);...
                  -2*normVec(1)*normVec(2), 1-2*normVec(2)^2];

    % Mode spectra, 1D curve, fine cubic interpolation
    Kx_list = [-0.1:0.002:0.1];
    Ky_list = [-0.1:0.002:0.1];
    [Kx_list1,Ky_list1] = meshgrid(Kx_list,Ky_list);
    Kx_list2 = reshape(Kx_list1,1,[]);
    Ky_list2 = reshape(Ky_list1,1,[]);

    Wpperr(Echeck==1)=0;
    Wmmerr(Echeck==1)=0;
    Wp_interp = griddata(K(:,1), K(:, 2), Wp2, Kx_list2, Ky_list2,'cubic');
    Wm_interp = griddata(K(:,1), K(:, 2), Wm2, Kx_list2, Ky_list2,'cubic');
    Wpe_interp = griddata(K(:,1), K(:, 2), sqrt(Wpperr),Kx_list2, Ky_list2,'cubic');
    Wme_interp = griddata(K(:,1), K(:, 2), sqrt(Wmmerr), Kx_list2, Ky_list2,'cubic');
    Wp_interp(Wp_interp<0)=0;
    Wm_interp(Wm_interp<0)=0;

    K2 = [Kx_list2',Ky_list2'];
    
    SegEig = zeros(size(SegDist,1),4);
    SegK_mirror = SegK*ReflectMat;
    for i=1:size(SegEig,1)
        % Consider rhombic symmetry, take average with mirrored K space.
        pos = [SegK(i,1),SegK(i,2)];
        pos_mirror = [SegK_mirror(i,1),SegK_mirror(i,2)];
        
        d = sqrt(sum((K2-pos).^2,2));
        d_mirror = sqrt(sum((K2-pos_mirror).^2,2));
        ind = find(d==min(d),1);
        ind_mirror = find(d_mirror==min(d_mirror),1);
        
        SegEig(i,1) = mean([Wp_interp(ind),Wp_interp(ind_mirror)]);
        SegEig(i,2) = mean([Wpe_interp(ind),Wpe_interp(ind_mirror)]);
        SegEig(i,3) = mean([Wm_interp(ind),Wm_interp(ind_mirror)]);
        SegEig(i,4) = mean([Wme_interp(ind),Wme_interp(ind_mirror)]);
    end
    
    % ANG model
    SegEigANG=zeros(size(SegDist));
    for i=1:size(SegK,1)
        SegEigANG_lower(i)=sqrt(lm(SegK(i,1),SegK(i,2)));
        SegEigANG_upper(i)=sqrt(lp(SegK(i,1),SegK(i,2)));
    end
    
    figure(4); %clf; 
    subplot(3,4,[10,11]);
    hold on; box on;
    %%%%%%%%%%%%%
    %  Frequency with errorbar
    errorbar(SegDist,SegEig(:,1),SegEig(:,2),'vertical','o')
    errorbar(SegDist,SegEig(:,3),SegEig(:,4),'vertical','o')
    %%%%%%%%%%%%%
    plot(SegDist,SegEigANG_lower,'-','LineWidth',2,'DisplayName','ANG model, lower')
    plot(SegDist,SegEigANG_upper,'-','LineWidth',2,'DisplayName','ANG model, upper')
    ylabel('Frequency'); ylim([-0.5e-3,4]); %axis auto
    currentaxis = axis();currentaxis(3)=0;
    axis(currentaxis);
    xticks(SegPtDist); xticklabels(SegPtName); xlim([0,max(SegPtDist)]);
    for i = 1:size(SegPtDist,1); xline(SegPtDist(i)); end
    legend
    legend({'Expt. mode lower','Expt. mode upper','ANG model lower','ANG model upper'},...
        'Location','southeast');
    
    subplot(3,4,12);
    hold off; 
    plot(SegPt(:,1),SegPt(:,2),'k-','LineWidth',2)
    hold on;
    plot(BoxX,BoxY,'k');
    axis equal; axis off;
    for i = [1,2,3,5]
        switch i
            case 1
                text(SegPt(i,1)-0.01,SegPt(i,2),SegPtName{i},'FontSize',15)
            otherwise
                text(SegPt(i,1)+0.005,SegPt(i,2),SegPtName{i},'FontSize',15)
        end
    end
    xlim([-0.07,0.07]); ylim([-0.07,0.07]);
        



