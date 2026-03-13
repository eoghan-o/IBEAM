function [qp, w] = tri_gauss3(T)
    A = polyarea(T(:,1), T(:,2));

    bary = [1/6 1/6 2/3;
        1/6 2/3 1/6;
        2/3 1/6 1/6];

    qp = bary * T;
    w  = (A/3)*ones(3,1);
end