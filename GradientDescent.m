function x = GradientDescent(f, x0, step_size, iterations)
%GRADIENTDESCENT Summary of this function goes here
%   Detailed explanation goes here
    x = x0(:);
    for i=1:iterations
        size(x)
        g = FiniteDiffGradient(f, x);
        x = x - step_size * g;
    end
end

