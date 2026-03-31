function B = bernoulli(p)
    if abs(p) < 1e-8
        B = 1 - p/2 + p^2;
    else
        B = p / (exp(p) - 1);
    end
end