function [qp, w] = tri_gauss9(T)
% TRI_GAUSS9 - 9-point Gaussian quadrature on a triangle
% 
% Inputs:
%   T : 3x2 array of triangle vertices [x y]
%       T = [x1 y1; x2 y2; x3 y3]
%
% Outputs:
%   qp : 9x2 array of quadrature points [x y]
%   w  : 9x1 vector of quadrature weights

% Area of the triangle
A = polyarea(T(:,1), T(:,2));

% Quadrature points in barycentric coordinates
alpha = 0.124949503233232;
beta  = 0.437525248383384;
gamma = 0.188203535619033;
delta = 0.810014365005957;

bary = [
    1/3,   1/3,   1/3;      % 1
    alpha, beta, 1-alpha-beta;  % 2
    beta,  alpha, 1-alpha-beta; % 3
    beta, 1-alpha-beta, alpha;  % 4
    gamma, gamma, 1-2*gamma;    % 5
    gamma, 1-2*gamma, gamma;    % 6
    1-2*gamma, gamma, gamma;    % 7
    delta, delta, 1-2*delta;    % 8
    delta, 1-2*delta, delta;    % 9
];

% Quadrature weights
w = A * [
    0.225;
    0.132394152788506;
    0.132394152788506;
    0.132394152788506;
    0.125939180544827;
    0.125939180544827;
    0.125939180544827;
    0.132394152788506;
    0.132394152788506;
];

% Convert barycentric to Cartesian coordinates
qp = bary * T;

end
