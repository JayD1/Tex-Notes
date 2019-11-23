function [x, R, P, Alpha, Beta] = cg(A,b,tol)
% CONJGRAD  Conjugate Gradient Method.
%   X = CONJGRAD(A,B) attemps to solve the system of linear equations A*X=B
%   for X. The N-by-N coefficient matrix A must be symmetric and the right
%   hand side column vector B must have length N.
%
%   X = CONJGRAD(A,B,TOL) specifies the tolerance of the method. The
%   default is 1e-10.
%
% Example (highlight lines between %{ and %}, then press F9):
%{
  n = 6000;
  m = 8000;
  A = randn(n,m);
  A = A * A';
  b = randn(n,1);
  tic, x = conjgrad(A,b); toc
  norm(A*x-b)
%}
% By Yi Cao at Cranfield University, 18 December 2008
% Updated on 6 Feb 2014.
%
    if nargin<3
        tol=1e-10;
    end
    x = b;
    r = b - A*x;
    if norm(r) < tol
        return
    end
%     r=r/norm(r);
    p = r;
    P = p;
    R = r;
    z = A*p;
    alpha = (r'*r)/(p'*z);
    Alpha = alpha;
    x = x + alpha*p;
    Beta = [0];
    for k = 1:numel(b);
       r = r - alpha*z;
       R = [R, r];
       beta = (r'*z)/(p'*z);
       Beta = [Beta, beta];
       p = r - beta*p;
       P = [P, p];
       z = A*p;
       alpha = (r'*p)/(p'*z);
       Alpha = [Alpha,alpha];
       x = x + alpha*p;
       if( norm(r) < tol )
            return;
       end
    end
end
 

