function [FEM, deleteBar] = Sensi(Phi, EW, EH, E, nu, Bars, RBF)
%=weipeng 2007-01-12=======================================================
% Finite element analysis
%-input parameters---------------------------------------------------------
% Phi       ...     the level set function value (nelx+1 by nely+1)
% EW        ...     the width of each element
% EH        ...     the height of each element
% idxWeak   ...     indicates if use weak material, 1: use; 0: don't use
% E         ...     
% nu        ...     
% Bars      ...     the coordinates of the bars [x1 y1 x2 y2 A1; x1 y1 x2 y2 A2;...]
% RBF       
%   .cent           the coordinates of the knots.
%   .cMQ            the c parameter for MQ RBF
%   .csWidth        the support width for CS RBF
%   .nKnots         the number of knots
%   .A              the collocation matrix
%   .invA           the inverse of the collocation matrix
% 	.coeff          the coefficients
%   .type           tpye of RBFs                                            % 1:IMQ, 2:C2 
%-output parameters--------------------------------------------------------
% FEM       ...
%   .com            the mean compliance of the structure
%   .den            the density of each element
%   .U              the displacement of nodes
%   .BSenX          the sensitivity of the objective function
%   .BConX          the sensitivity of the constraint
%   .BVolume        the volume of bars
% deleteBar ...     the tag for the bars should be deleted in the next loop
%==========================================================================
deleteBar = []; 
[nely0,nelx0]=size(Phi);
nely=nely0-1;
nelx=nelx0-1;

h = 1;
[KB] = BasicKe(E,nu, EW, EH,h);
den = ones(nely, nelx);

K = sparse(2*(nelx+1)*(nely+1), 2*(nelx+1)*(nely+1));
F = sparse(2*(nely+1)*(nelx+1),1);
U = zeros(2*(nely+1)*(nelx+1),1);
n = 0;
for elx = 1:nelx
    for ely = 1:nely
        n=n+1;
        n1 = (nely+1)*(elx-1)+ely;
        n2 = (nely+1)* elx   +ely;
        edof = [2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n2+1; 2*n2+2; 2*n1+1; 2*n1+2];
      	K(edof,edof) = K(edof,edof) + den(ely,elx)*KB;
    end
end
                                                                          
% === codimension problems with bar ===
BarE = E*1000;
[nBar,~]= size(Bars);
BarXY = Bars(:,1:4);                                                        % the coordinates for bars in 2 dimensional
BarA = Bars(:,5);
BSenX = zeros(nBar,5);                                                      % sensitivity of obj
BConX = zeros(nBar,5);                                                      % sensitivity of constraints

% === for RBF ===
KeBarN = zeros(2*(nely+1)*(nelx+1),2*(nely+1)*(nelx+1),nBar);

T1 = zeros(4,2*(nelx+1)*(nely+1));
T2 = zeros(2*(nelx+1)*(nely+1),2*(nelx+1)*(nely+1));
T2(1:(nelx+1)*(nely+1),1:(nelx+1)*(nely+1))         = RBF.invA;
T2((nelx+1)*(nely+1)+1:end,(nelx+1)*(nely+1)+1:end) = RBF.invA;

for iBar = 1:nBar
    DX = BarXY(iBar,1)-BarXY(iBar,3);
    DY = BarXY(iBar,2)-BarXY(iBar,4);
    lBar(iBar) = sqrt(DX^2+DY^2);                                           % the length of the bar

    Ke11 = [DX*DX DX*DY; DX*DY DY*DY];
    pKeBarpA = BarE/(lBar(iBar)^3)*[Ke11 -Ke11; -Ke11 Ke11];                % for sensitivity analysis of area (special)
    KeBar(:,:,iBar) = BarA(iBar)*pKeBarpA;

    % === for RBF ===
    edof1 = (1 : (nelx+1)*(nely+1));
    edof2 = ((nelx+1)*(nely+1)+1 :2*(nelx+1)*(nely+1));
    Kedof = [edof1;edof2];
    Kedof = Kedof(:);

    % == RBF ===
    TPhiNode1 = DirectMatrix(RBF.cent,[BarXY(iBar,1),BarXY(iBar,2)],RBF.type,RBF.cMQ,RBF.csWidth);
    TPhiNode2 = DirectMatrix(RBF.cent,[BarXY(iBar,3),BarXY(iBar,4)],RBF.type,RBF.cMQ,RBF.csWidth);
    T1 = T1*0;
    T1(1,1:(nelx+1)*(nely+1))       = TPhiNode1;
    T1(2,(nelx+1)*(nely+1)+1:end)   = TPhiNode1;
    T1(3,1:(nelx+1)*(nely+1))       = TPhiNode2;
    T1(4,(nelx+1)*(nely+1)+1:end)   = TPhiNode2;
    T = T1*T2;

    KeBarN(:,:,iBar) = T(:,:)'*KeBar(:,:,iBar)*T(:,:);
    K = K + KeBarN(Kedof(:), Kedof(:),iBar);                                % for RBF
end

% DEFINE LOADS AND SUPPORTS (Cantilever-BEAM)
ny=ceil(nely/2);
F(2*((nely+1)*nelx+(ny+1)),1) = -1e5;
fixeddofs   = (1:1:2*(nely+1));
alldofs     = (1:2*(nely+1)*(nelx+1));
freedofs    = setdiff(alldofs,fixeddofs);

U(freedofs,:) = K(freedofs,freedofs) \ F(freedofs,:);
U(fixeddofs,:)= 0;
% === for sensitivity analysis === weipeng: 2008-09-08
for iBar = 1:nBar
    DeltaX1 = BarXY(iBar,1)-BarXY(iBar,3);
    DeltaY1 = BarXY(iBar,2)-BarXY(iBar,4);
    DeltaX2 = BarXY(iBar,3)-BarXY(iBar,1);
    DeltaY2 = BarXY(iBar,4)-BarXY(iBar,2);
    BUe   = U(:,1);

    Ke11 = [DeltaX1*DeltaX1 DeltaX1*DeltaY1; DeltaX1*DeltaY1 DeltaY1*DeltaY1];
    pKeBarpA = BarE/(lBar(iBar)^3)*[Ke11 -Ke11; -Ke11 Ke11];                % for sensitivity analysis of area (special)

    SampleNode = [BarXY(iBar,1) BarXY(iBar,2);BarXY(iBar,3) BarXY(iBar,4)];
    [Phix, pPhipx, pPhipy]=PartialDirectMatrix(RBF.cent, SampleNode, RBF.type, RBF.cMQ, RBF.csWidth);

    TPhiNode1 = Phix(1,:);
    TPhiNode2 = Phix(2,:);

    T1 = T1*0;
    T1(1,1:(nelx+1)*(nely+1))       = TPhiNode1;
    T1(2,(nelx+1)*(nely+1)+1:end)   = TPhiNode1;
    T1(3,1:(nelx+1)*(nely+1))       = TPhiNode2;
    T1(4,(nelx+1)*(nely+1)+1:end)   = TPhiNode2;
%     T = T1*T2;
    
    pTpx1 = T1*0;                                                           % partial derivative respect to x1
    pTpx1(1,1:(nelx+1)*(nely+1))       = pPhipx(1,:);
    pTpx1(2,(nelx+1)*(nely+1)+1:end)   = pPhipx(1,:);

    pTpy1 = T1*0;                                                           % partial derivative respect to y1
    pTpy1(1,1:(nelx+1)*(nely+1))       = pPhipy(1,:);
    pTpy1(2,(nelx+1)*(nely+1)+1:end)   = pPhipy(1,:);

    pTpx2 = T1*0;                                                           % partial derivative respect to x2
    pTpx2(3,1:(nelx+1)*(nely+1))       = pPhipx(2,:);
    pTpx2(4,(nelx+1)*(nely+1)+1:end)   = pPhipx(2,:);

    pTpy2 = T1*0;                                                           % partial derivative respect to y2
    pTpy2(3,1:(nelx+1)*(nely+1))       = pPhipy(2,:);
    pTpy2(4,(nelx+1)*(nely+1)+1:end)   = pPhipy(2,:);

    % === partial KE partial X 1234 ===
    Kx1 = [2*DeltaX1 DeltaY1; DeltaY1 0         ];    Kx2 = [2*DeltaX2 DeltaY2; DeltaY2 0           ];
    Ky1 = [0         DeltaX1; DeltaX1 2*DeltaY1 ];    Ky2 = [0         DeltaX2; DeltaX2 2*DeltaY2   ];
    pKpx1 = -3*KeBar(:,:,iBar)*DeltaX1/(lBar(iBar)^2) + BarE*BarA(iBar)/(lBar(iBar)^3)*[Kx1 -Kx1; -Kx1 Kx1];
    pKpy1 = -3*KeBar(:,:,iBar)*DeltaY1/(lBar(iBar)^2) + BarE*BarA(iBar)/(lBar(iBar)^3)*[Ky1 -Ky1; -Ky1 Ky1];
    pKpx2 = -3*KeBar(:,:,iBar)*DeltaX2/(lBar(iBar)^2) + BarE*BarA(iBar)/(lBar(iBar)^3)*[Kx2 -Kx2; -Kx2 Kx2];
    pKpy2 = -3*KeBar(:,:,iBar)*DeltaY2/(lBar(iBar)^2) + BarE*BarA(iBar)/(lBar(iBar)^3)*[Ky2 -Ky2; -Ky2 Ky2];
    % KeBarN(:,:,iBar) = T(:,:)'*KeBar(:,:,iBar)*T(:,:);

    pKBpx1(:,:) = 2*pTpx1'*KeBar(:,:,iBar)*T1 + T1'*pKpx1*T1;
    pKBpy1(:,:) = 2*pTpy1'*KeBar(:,:,iBar)*T1 + T1'*pKpy1*T1;
    pKBpx2(:,:) = 2*pTpx2'*KeBar(:,:,iBar)*T1 + T1'*pKpx2*T1;
    pKBpy2(:,:) = 2*pTpy2'*KeBar(:,:,iBar)*T1 + T1'*pKpy2*T1;
    
    BUe = T2(Kedof,Kedof)*BUe;
    
    BSenX(iBar,1) = -1/2 * (BUe' * pKBpx1(Kedof,Kedof) * BUe);
    BSenX(iBar,2) = -1/2 * (BUe' * pKBpy1(Kedof,Kedof) * BUe);
    BSenX(iBar,3) = -1/2 * (BUe' * pKBpx2(Kedof,Kedof) * BUe);
    BSenX(iBar,4) = -1/2 * (BUe' * pKBpy2(Kedof,Kedof) * BUe);

    pKeBarpA = T1'*pKeBarpA*T1;
    BSenX(iBar,5) = -1/2 * (BUe' * pKeBarpA(Kedof,Kedof) * BUe);            % sensitivity of area (special)

    BConX(iBar,1) = Bars(iBar,5)*DeltaX1/lBar(iBar);                        % Volume derivate with respect to X1
    BConX(iBar,2) = Bars(iBar,5)*DeltaY1/lBar(iBar);
    BConX(iBar,3) = Bars(iBar,5)*DeltaX2/lBar(iBar);
    BConX(iBar,4) = Bars(iBar,5)*DeltaY2/lBar(iBar);                        % Volume derivate with respect to Y2
    BConX(iBar,5) = lBar(iBar);                                             % Volume derivate with respect to Area
end

c = 1/2*U'*K*U;
% === deleting bars ===
BVolume = 0;
for iBar = 1:nBar
    iVolume = lBar(iBar)*Bars(iBar,5);
    if lBar(iBar) <= 1e-4 || Bars(iBar,5) <= 0.01
        deleteBar = [deleteBar;iBar];
    end
    BVolume = BVolume + iVolume;
end 
FEM.com = c;
FEM.den = den;
FEM.U = U;
FEM.BSenX = BSenX;
FEM.BConX = BConX;
FEM.BVolume = BVolume;

function [KE] = BasicKe(E,nu, a, b,h)
k = [-1/6/a/b*(nu*a^2-2*b^2-a^2),  1/8*nu+1/8, -1/12/a/b*(nu*a^2+4*b^2-a^2), 3/8*nu-1/8, ...
      1/12/a/b*(nu*a^2-2*b^2-a^2),-1/8*nu-1/8,  1/6/a/b*(nu*a^2+b^2-a^2),   -3/8*nu+1/8];
KE = h*E/(1-nu^2)*...
                [ k(1) k(2) k(3) k(4) k(5) k(6) k(7) k(8)
                  k(2) k(1) k(8) k(7) k(6) k(5) k(4) k(3)
                  k(3) k(8) k(1) k(6) k(7) k(4) k(5) k(2)
                  k(4) k(7) k(6) k(1) k(8) k(3) k(2) k(5)
                  k(5) k(6) k(7) k(8) k(1) k(2) k(3) k(4)
                  k(6) k(5) k(4) k(3) k(2) k(1) k(8) k(7)
                  k(7) k(4) k(5) k(2) k(3) k(8) k(1) k(6)
                  k(8) k(3) k(2) k(5) k(4) k(7) k(6) k(1)];        

