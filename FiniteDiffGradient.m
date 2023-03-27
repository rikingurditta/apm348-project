function g = FiniteDiffGradient(f, x)
%FINITEDIFFGRADIENT Summary of this function goes here
%   Detailed explanation goes here
    h = 0.00001;
    n = size(x);
    g = zeros(n);
    for i=1:n
        x1 = x(:);
        x2 = x(:);
        x1(i) = x(i) + h;
        x2(i) = x(i) - h;
        g(i) = (f(x1) - f(x2)) / h;
    end
end

