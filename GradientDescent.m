function x = GradientDescent(f, x0, step_size, iterations)
%GRADIENTDESCENT use FiniteDiffGradient to perform gradient descent on f
    % start at initial guess
    x = x0(:);
    for i=1:iterations
        display(['iteration ', num2str(i)])
        display(['x ', num2str(x')])
        display(['f(x) ', num2str(f(x))])
        g = FiniteDiffGradient(f, x);
        display(['g ', num2str(g')])
        fprintf('\n');
        x = x - step_size * g;
    end
end

