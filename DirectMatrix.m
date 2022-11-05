function [AA]=DirectMatrix(knot,sampleCoord,rType,cMQ,csWidth)
[nColumn,dim]=size(knot);
[nLine,dim]=size(sampleCoord);
A=sparse(nLine,nColumn); 
AA=zeros(nLine,nColumn);
for i=1:dim
    A=A+( ones(nLine,1)*knot(:,i)' - sampleCoord(:,i)*ones(1,nColumn) ).^2;
end
c = csWidth;
switch rType
    case 1  %IMQ
        A = A + cMQ^2;     %now, RBF.c can only be a scalar, special dealing should be done to make it a matrix.
        AA = A.^(-0.5);

    case 2  %CSRBF2 
        A=(A.^(0.5))./ c;
        ind=find(A<1);  % indices of entries which are less than zero.
        AA(ind)=(( 1- A(ind) ).^4) .* ( 4*A(ind) + 1 );

    otherwise
        error('Wrong type of RBF!');
end
