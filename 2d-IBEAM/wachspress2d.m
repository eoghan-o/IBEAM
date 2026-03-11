function [phi, dphi] = wachspress2d(v,x)
% Written by M. Floater, et al.
% Evaluate Wachspress basis functions and their gradients in a convex
% polygon
%
% Inputs:
% v     : [x1 y1; x2 y2; ...; xn yn], the n vertices of the polygon in ccw
% x     : [x(1) x(2)], the point at which the basis functions are computed
% Outputs:
% phi   : output basis functions = [phi_1; ...; phi_n]
% dphi  : output gradient of basis functions = [dphi_1; ...; dphi_n]
    n = size(v,1);
    w = zeros(n,1);
    R = zeros(n,2);
    phi = zeros(n,1);
    dphi = zeros(n,2);

    un = getNormals(v);
    p = zeros(n,2);
    for i = 1:n
        h = dot(v(i,:) - x,un(i,:));
        p(i,:) = un(i,:) / h;
    end

    for i = 1:n
        im1 = mod(i-2,n) + 1;
        w(i) = det([p(im1,:); p(i,:)]);
        R(i,:) = p(im1,:) + p(i,:);
    end

    wsum = sum(w);
    phi = w/wsum;

    phiR = phi' * R;
    for k = 1:2
        dphi(:,k) = phi .* (R(:,k) - phiR(:,k));
    end
end

function un = getNormals(v)
% Written by M. Floater, et al.
% Function to compute the outword unit normal to each edge
    n = size(v,1);
    un = zeros(n,2);
    for i = 1:n
        d = v(mod(i,n) + 1,:) - v(i,:);
        un(i,:) = [d(2) -d(1)]/norm(d);
    end
end