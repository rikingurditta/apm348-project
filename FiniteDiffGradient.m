function g = FiniteDiffGradient(f, x)
%FINITEDIFFGRADIENT use finite differences to find gradient of f
    h = 0.000000001;
    n = size(x);
    g = zeros(n);
    % for each i, get i'th partial derivative
    for i=1:n
        x1 = x(:);
        x2 = x(:);
        x1(i) = x(i) + h;
        x2(i) = x(i) - h;
        % central finite difference
        g(i) = (f(x1) - f(x2)) / (2 * h);
    end
end

