function x_hat = AMP(y, A, alhpa)
    [M,N] = size(A);
    assert(M == length(y) && isvector(y));
    A_H = A';
    x_hat = zeros(N,1); 
    v = zeros(M,1); d = 0; dc = 0;
    t = 0;
    while t <= 40
        v = y-A*x_hat+(d/M)*v+(dc/M)*conj(v);
        r = x_hat+A_H*v;
        sigma2 = norm(v)^2/M;
        [x_hat, d, dc] = shrink(r, alhpa, sigma2);
        t = t+1;
    end 
end
function [x_hat, d, dc] = shrink(r, alpha, sigma2)
    sigma = sqrt(sigma2);
    b = abs(r)>alpha*sigma;
    x_hat = (r-alpha*sigma*exp(1j*angle(r))).*b;
    d = sum((1-alpha*sigma./(2*abs(r))).*b);
    dc = -sum(alpha*sigma*(r.^2)./(2*abs(r).^3).*b);
end

