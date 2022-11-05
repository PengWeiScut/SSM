function [AA, pApx, pApy]=PartialDirectMatrix(knot,sampleCoord,rType,cMQ,csWidth)

[nColumn,~]=size(knot);
[nLine,dim]=size(sampleCoord);
A=sparse(nLine,nColumn); 
AA=zeros(nLine,nColumn);
pA=zeros(nLine,nColumn);
pApx=zeros(nLine,nColumn);
pApy=zeros(nLine,nColumn);
for i=1:dim
    A=A+( ones(nLine,1)*knot(:,i)' - sampleCoord(:,i)*ones(1,nColumn) ).^2;
end
c = csWidth;  % matrix of absolute shape parameter
switch rType
    case 1  %IMQ
        A = A + cMQ^2;     %now, RBF.c can only be a scalar, special dealing should be done to make it a matrix.
        AA = A.^(-0.5);
        pApx = ( ones(nLine,1)*knot(:,1)' - sampleCoord(:,1)*ones(1,nColumn) ) .*(AA.^3);
        pApy = ( ones(nLine,1)*knot(:,2)' - sampleCoord(:,2)*ones(1,nColumn) ) .*(AA.^3);
    case 2  %CSRBF2
        A=(A.^(0.5))./ c;
        ind=find(A<1);     % indices of entries which are less than zero.
        pA(ind)=(( 1- A(ind) ).^3) .* ( -20*A(ind) );
        pApx = - pA .* ( ones(nLine,1)*knot(:,1)' - sampleCoord(:,1)*ones(1,nColumn) ) ./ (A + realmin)/c/c;
        pApy = - pA .* ( ones(nLine,1)*knot(:,2)' - sampleCoord(:,2)*ones(1,nColumn) ) ./ (A + realmin)/c/c;
        AA(ind)=(( 1- A(ind) ).^4) .* ( 4*A(ind) + 1 );
    otherwise
        error('Wrong type of RBF!');
end
