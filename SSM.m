function []=SSM()
%==========================================================================
% The program, SSM, is developed for the implimentation of  topology
% optimization of truss structures with Stiffness Spreading Method. 
% First Version Date: Apr. 30, 2010 
% Current Version: 2022.11.04.001
%
% Developed by: Wei Peng
% Email: ctpwei@scut.edu.cn
%
% Main references:
% Peng Wei, Haitao Ma, Michael Yu Wang, The Stiffness Spreading Method 
% for Layout Optimization of Truss Structures, Structural and 
% Multidisciplinary Optimization, 49(2), 667-682, April 2014. 
%
% Peng Wei, Haitao Ma, Taicong Chen and Michael Yu Wang, Stiffness 
% Spreading Method for Layout Optimization of Truss Structures, 
% 6th China-Japan-Korea Joint Symposium on Optimization of Structural
% and Mechanical Systems, June 22¨C25, 2010, Kyoto, Japan 
% 
% **************************   Disclaimer   ***************************** %
% The  authors  reserve all rights for the programs. The programs  may  be
% distributed  and used for academic and educational  purposes. The authors
% do not guarantee that the code is free from errors, and they shall not be
% liable in any event caused by the use of the programs.
%==========================================================================

%% Fundamental parameters
Domain.nelx = 20;                                        
Domain.nely = 10;
nelx = Domain.nelx;
nely = Domain.nely;
Domain.size = [nelx nely];                              
Domain.size(2) = Domain.size(1)*nely/nelx;              
Domain.aspect = Domain.size(1)/nelx;                     
EleWidth = Domain.aspect; EleHeight = Domain.aspect;    
ti1 = 0 : EleWidth : Domain.size(1);
ti2 = 0 : EleHeight : Domain.size(2);
[XI,YI] = meshgrid(ti1,ti2);                             
xD = XI(:);
yD = YI(:);
Domain.nodeGrid=[xD yD];
Phi = ones(Domain.nelx+1, Domain.nely+1);              
E = 2.1e11; NU = 0.3;                                                                          
% === Initialize the initial design === %
Bar1 = [ 0 0 nelx/10 nely/5;
         0 nely/5 nelx/10 0 ];                                              % one cross
Bars = ones(20,5)*1;                                 
for i = 1:3  
    Barx = Bar1;
    Barx(:,[1 3]) = Bar1(:,[1 3]) + 0.2*nelx + (i-1)*nelx*0.25;
    Barx(:,[2 4]) = Bar1(:,[2 4]) + 0.1*nely;                
    Bars(i*2-1:i*2,1:4)=Barx;                         
end
for i = 4:7  
    j = i-3;                                       
    Barx = Bar1;
    Barx(:,[1 3]) = Bar1(:,[1 3]) + 0.075*nelx + (j-1)*nelx*0.25;
    Barx(:,[2 4]) = Bar1(:,[2 4]) + 0.4*nely;
    Bars(i*2-1:i*2,1:4)=Barx;
end
for i = 8:10 
    j = i-7;
    Barx = Bar1;
    Barx(:,[1 3]) = Bar1(:,[1 3]) + 0.2*nelx + (j-1)*nelx*0.25;      
    Barx(:,[2 4]) = Bar1(:,[2 4]) + 0.7*nely;
    Bars(i*2-1:i*2,1:4)=Barx;
end
[nBar,~]= size(Bars);                                                                        
BVolume = 0;                                              
for iBar = 1:nBar
    DX = Bars(iBar,1)-Bars(iBar,3);
    DY = Bars(iBar,2)-Bars(iBar,4);
    lBar = sqrt(DX^2+DY^2);                                                 % the length of the bar
    BVolume = BVolume + lBar*Bars(iBar,5);                
end

%% RBF Initialization
RBF.cent=Domain.nodeGrid;
RBF.cMQ = 3;
RBF.csWidth = 2.2;
RBF.nKnots = length(Domain.nodeGrid);                     
RBF.A = 0;
RBF.invA = 0;
RBF.coeff = 0;
RBF.type = 2;                                                               % 1:IMQ, 2:C2 
% === initial RBF coefficients === %
np = (nelx+1)*(nely+1);                                                     % number of points
if 2 == RBF.type 
	RBF.A = sparse(np,np);                                 
end
[RBF.A]=DirectMatrix(RBF.cent,RBF.cent,RBF.type,RBF.cMQ ,RBF.csWidth);
RBF.invA = inv(RBF.A); 
%% The parameters used for MMA
m = 1;                                                                      % number of the contraint
n = nBar*5;                                                                 % number of the design variables
Amax = 3;
xmin  = 0.0001*ones(n,1);
xmax  = ones(nBar,1);
xmax = kron([Domain.size Domain.size],xmax)';
xmax = [xmax; ones(nBar,1)'*Amax];
xmax = xmax(:); 
xval = Bars(:,1:5)';                                                        % Alpha, and the design variable
xval = xval(:);                                      
low   = xmin; upp   = xmax; 

c = 1e6*ones(m,1);                      d = ones(m,1);
a0 = 1;                                 
a = 1e-3*ones(m,1);
%% Optimization Initialization
OPT.m = m;                                                                  % the number of constraint
OPT.n = n;                                                                  % the number of design variable

OPT.xval  = xval;
OPT.f0val = 0;                                                              % Objective function
OPT.fval = 0;                                                               % constraint function
OPT.xold1 = xval;
OPT.xold2 = xval;
OPT.dfdx   = zeros(n,1);
OPT.dfdx2  = zeros(1,n);
OPT.df0dx  = zeros(n,1);
OPT.df0dx2 = zeros(n,1);
OPT.xmin = xmin;
OPT.xmax = xmax;
OPT.low = low;
OPT.upp = upp;

% === output window configuration === %
xAxis = 0:Domain.aspect:Domain.size(1);
yAxis = 0:Domain.aspect:Domain.size(2);
FigNum = 0;
set(0,'Units','pixels');
get(0,'ScreenSize');

%% Optimization Iteration
ItNum = 0;
TotalItNum = 500; 
change = 1;
while ItNum < TotalItNum && change > 1e-8
    ItNum = ItNum + 1;
    FigNum = FigNum + 1;
    disp(' ');
    disp(['FEM analysis No. ', num2str(ItNum), ' starts. ']);
    
    [FEM,deleteBar]=Sensi(Phi, Domain.aspect, Domain.aspect, E, NU, Bars, RBF); 
	FEM.BConX(:,1:4) = 0;
    % === delete bars === %
    nDelBar = length(deleteBar);
    if nDelBar > 0
        for iNum = 1:nDelBar
            iDelBar = deleteBar(nDelBar-iNum+1);                            % backward sequence is needed because the elements in matrix can only be deleted from the end to the beginning
            idx = (iDelBar-1)*5+1:iDelBar*5;
   
            OPT.xval(idx)  = [];
            OPT.xold1(idx)  = [];
            OPT.xold2(idx)  = [];
            OPT.dfdx2 (idx)  = [];
            OPT.df0dx2(idx)  = [];
            OPT.xmin(idx)  = [];
            OPT.xmax(idx)  = [];
            OPT.low(idx)  = [];
            OPT.upp(idx)  = [];
            FEM.BSenX(iDelBar,:) = [];
            FEM.BConX(iDelBar,:)  = [];    
        end
        Bars(deleteBar,:)  = [];
        n = n - nDelBar*5;
        OPT.n = n;
    end

    VolumeRatio(ItNum) = FEM.BVolume;
    Obj(ItNum) = FEM.com;
    Phi = reshape(Phi(:),nely+1, nelx+1);

%% Display 
    figure(1);    clf;
    hold on;
    surf(xAxis,yAxis,Phi*0,'FaceColor',[0.8 0.8 0.9]);

    % === Draw Bars === %
    [BNum,~]= size(Bars);
    line([0,0,nelx,nelx,0] , [0,nely,nely,0,0],'color',[0.3 0.2 0.8],'LineWidth',2);hold on;
    for i = 1:BNum
        line([Bars(i,1),Bars(i,3)],[Bars(i,2),Bars(i,4)],'color','r','LineWidth',Bars(i,5)*10);
        h_Nodes(i) = plot([Bars(i,1),Bars(i,3)],[Bars(i,2),Bars(i,4)],'o');
        set(h_Nodes(i),'LineWidth',1,'MarkerEdgeColor','w','MarkerFaceColor',[0.1 0.2 0.3]);
    end
    axis equal; axis off;
    
    figure(2);
    subplot(2,1,1);
    plot(log(Obj),'-rs','LineWidth',1,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','g',...
        'MarkerSize',1);
    xlabel('Iteration number');
    ylabel('ln(J)');
    title(['Total Number of Iterations: ',num2str(ItNum)]);
    axis([0 ItNum 0 1]);  axis 'auto y';

    subplot(2,1,2);
    plot(VolumeRatio,'-b*','LineWidth',1,...
        'MarkerEdgeColor','g',...
        'MarkerFaceColor','r',...
        'MarkerSize',1);
    xlabel('Iteration number');
    ylabel('Volume');
    axis([0 ItNum 0 1]);  axis 'auto y';
    pause(1e-6);

%% Bars update 
    [nBar,~]= size(Bars);
    tempBSenX = FEM.BSenX;
    tempBConX = FEM.BConX;
    OPT.f0val = Obj(ItNum);
    tempBSenX = tempBSenX';
    OPT.df0dx = tempBSenX(:);
    tempBConX = tempBConX';
    OPT.dfdx = tempBConX(:);
    OPT.fval = FEM.BVolume-BVolume;
    [xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,OPT.low,OPT.upp] = mmasub(m,n,ItNum,OPT.xval,OPT.xmin, OPT.xmax, ...
        OPT.xold1,OPT.xold2,OPT.f0val,OPT.df0dx,OPT.df0dx2,OPT.fval,OPT.dfdx,OPT.dfdx2,OPT.low,OPT.upp,a0,a,c,d);
    OPT.xold2 = OPT.xold1;     OPT.xold1 = OPT.xval;    OPT.xval = xmma;
    Bars(:,1:5) = reshape(OPT.xval,5,nBar)';
    if ItNum > 1 
       change = abs(Obj(ItNum)-Obj(ItNum-1));
    end
end
