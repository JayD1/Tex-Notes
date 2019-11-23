n=1000                                                                                                                                     ;
[U,S,V]=svd(randn(n));
s=diag(S);
A=U*diag(s+max(s))*U';
b=randn(n,1);
tic,[x,R,P,Alpha,Beta]=cg(A,b);toc
tic,x1=A\b;toc
norm(x-x1)
err = @(x) sqrt(sum(sum(x^2)));
norm(A*x-b)
[var,len] = size(P);
x = P'*A*P;
y = R'*R;
z = R'*A*R;
Z = z;
T = eye(len);
L = eye(len);
D = eye(len);
Lambda = [0];
Eta = [Alpha(1)];
for i=2:len
    Lambda = [Lambda, Beta(i)/Eta(i-1)];
    Eta = [Eta, Alpha(i)-Beta(i)*Lambda(i)];
end
for i=1:len
    x(i,i)=0;
    y(i,i)=0;
    z(i,i)=0;
    T(i,i) = Alpha(i);
    D(i,i) = Eta(i);
    if i ~= len
        L(i+1,i) = Lambda(i+1);
        z(i,i+1)=0;
        T(i,i+1) = Beta(i+1);
        z(i+1,i)=0;
        T(i+1,i) = Beta(i+1);
    end
end
err(x)
err(y)
err(z)