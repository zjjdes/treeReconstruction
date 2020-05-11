function s = svfun(d)
% 50 check points
    a1 = 1.182e+04;
    b1 = 0.8612;
    c1 = 1.004;
    range = 0.9;
    if d <= range
        s = a1*exp(-((d-b1)/c1)^2);
    else
        s = a1*exp(-((range-b1)/c1)^2);
    end
end