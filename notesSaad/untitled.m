n = 2000;
Q = RandOrthMat(n);
Lambda = diag(randi(n,n,1));
A = Q*Lambda*Q';
x=randi(n,n,1);
b=randi(n,n,1);
r = b-A*x;
p = r;
error = 0;
x_exact = inv(A)*b;
i=0;

while 1
    p_temp = A*p;
    alpha = (r'*r)/(p_temp'*p);
    x = x + alpha*p;
    temp = r'*r;
    r = r - alpha*p_temp;
    beta = (r'*r)/temp;
    p = r + beta*p;
    if sqrt(mean((x - x_exact).^2)) == error
        break
    end
    i=i+1
    error = sqrt(mean((x - x_exact).^2))
end

% error = sqrt(mean((x - x_exact).^2))