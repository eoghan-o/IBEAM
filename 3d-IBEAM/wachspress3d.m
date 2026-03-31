function [phi, dphi] = wachspress3d(v,g,un,x)
% Evaluate Wachspress basis functions and their gradients
% in a convex polyhedron
%
% Inputs:
% v : [x1 y1 z1; x2 y2 z2; ...; xn yn zn], the n vertices of the polyhedron
% g : cell array: [i1 i2 ... i_{k1}]; . . . ; [i1 i2 ... i_{kn}],
% which are the n neighborhood graphs, each in some counter-clockwise
% order as seen from outside the polyhedron
% un : [x1 y1 z1; x2 y2 z2; ...; xm ym zm], unit normal to each facet
% x : [x(1) x(2) x(3)], the point at which the basis functions are computed
% Outputs:
% phi : basis functions = [phi_1; ...; phi_n]
% dphi : gradient of basis functions = [dphi_1; ...; dphi_n]
    n = size(v,1);
    w = zeros(n,1);
    R = zeros(n,3);
    phi = zeros(n,1);
    dphi = zeros(n,3);
    for i = 1:n
        f = g{i};
        k = length(f);
        p = zeros(k,3);
        for j = 1:k
            h = dot(v(i,:) - x,un(f(j),:));
            p(j,:) = un(f(j),:) / h;
        end
        wloc = zeros(k-2,1);
        Rloc = zeros(k-2,3);
        for j = 1:k-2
            wloc(j) = det([p(j,:); p(j+1,:); p(k,:)]);
            Rloc(j,:) = p(j,:) + p(j+1,:) + p(k,:);
        end

        w(i) = sum(wloc);
        R(i,:) = (wloc' * Rloc) / w(i);
    end

    wsum = sum(w);
    phi = w/wsum;
    phiR = phi' * R;
    for d = 1:3
        dphi(:,d) = phi .* (R(:,d) - phiR(:,d));
    end
end