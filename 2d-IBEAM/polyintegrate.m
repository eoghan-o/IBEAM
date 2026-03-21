function total = polyintegrate(v,n)
    n_v = size(v, 1);
    total = 0;
    for i = 1:(n_v-2)
        j = i + 1;
        k = i + 2;
        t = v([1 j k], :);
        total = total + triquad(t,n);
    end
end