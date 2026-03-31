%% test_polyintegrate3d.m
% Tests polyintegrate3d against known analytical integrals
% on a unit cube and a reference tetrahedron

fprintf('=== Testing polyintegrate3d ===\n\n');

%% -------------------------------------------------------
% TEST 1: Unit tetrahedron, integrand = 1
% Exact volume = 1/6
% -------------------------------------------------------
v_tet = [0 0 0;
         1 0 0;
         0 1 0;
         0 0 1];

f_tet = {[2 3 4],   % base (does not contain vertex 1)
         [1 3 4],   % skip (contains vertex 1)
         [1 2 4],   % skip
         [1 2 3]};  % skip

n_quad = 4;
result1 = polyintegrate3d(v_tet, f_tet, n_quad, @(x) 1);
exact1  = 1/6;
fprintf('Test 1 - Unit tet, f=1 (volume)\n');
fprintf('  Result: %.10f\n', result1);
fprintf('  Exact:  %.10f\n', exact1);
fprintf('  Error:  %.2e\n\n', abs(result1 - exact1));

%% -------------------------------------------------------
% TEST 2: Unit tetrahedron, integrand = x
% Exact = integral of x over unit tet = 1/24
% -------------------------------------------------------
result2 = polyintegrate3d(v_tet, f_tet, n_quad, @(x) x(1));
exact2  = 1/24;
fprintf('Test 2 - Unit tet, f=x\n');
fprintf('  Result: %.10f\n', result2);
fprintf('  Exact:  %.10f\n', exact2);
fprintf('  Error:  %.2e\n\n', abs(result2 - exact2));

%% -------------------------------------------------------
% TEST 3: Unit tetrahedron, integrand = x*y*z
% Exact = 1/720
% -------------------------------------------------------
result3 = polyintegrate3d(v_tet, f_tet, n_quad, @(x) x(1)*x(2)*x(3));
exact3  = 1/720;
fprintf('Test 3 - Unit tet, f=xyz\n');
fprintf('  Result: %.10f\n', result3);
fprintf('  Exact:  %.10f\n', exact3);
fprintf('  Error:  %.2e\n\n', abs(result3 - exact3));

%% -------------------------------------------------------
% TEST 4: Unit cube [0,1]^3, integrand = 1
% Exact volume = 1
% Cube split into 6 tetrahedra (standard decomposition)
% -------------------------------------------------------
v_cube = [0 0 0;   % 1
          1 0 0;   % 2
          1 1 0;   % 3
          0 1 0;   % 4
          0 0 1;   % 5
          1 0 1;   % 6
          1 1 1;   % 7
          0 1 1];  % 8

% 6-tet decomposition of unit cube
% Each row is [v1 v2 v3 v4] of one tet
tets = [1 2 4 5;
        2 4 5 8;
        2 5 6 8;
        2 3 4 8;
        2 6 7 8;
        2 3 7 8];

result4 = zeros(1,1);
for t = 1:6
    vt = v_cube(tets(t,:), :);
    ft = {[2 3 4],   % only face not containing local vertex 1
          [1 3 4],
          [1 2 4],
          [1 2 3]};
    result4 = result4 + polyintegrate3d(vt, ft, n_quad, @(x) 1);
end
exact4 = 1.0;
fprintf('Test 4 - Unit cube, f=1 (volume)\n');
fprintf('  Result: %.10f\n', result4);
fprintf('  Exact:  %.10f\n', exact4);
fprintf('  Error:  %.2e\n\n', abs(result4 - exact4));

%% -------------------------------------------------------
% TEST 5: Unit cube, integrand = x^2 + y^2 + z^2
% Exact = integral over [0,1]^3 of x^2+y^2+z^2 = 3*(1/3) = 1
% -------------------------------------------------------
result5 = zeros(1,1);
for t = 1:6
    vt = v_cube(tets(t,:), :);
    ft = {[2 3 4],
          [1 3 4],
          [1 2 4],
          [1 2 3]};
    result5 = result5 + polyintegrate3d(vt, ft, n_quad, @(x) x(1)^2 + x(2)^2 + x(3)^2);
end
exact5 = 1.0;
fprintf('Test 5 - Unit cube, f=x^2+y^2+z^2\n');
fprintf('  Result: %.10f\n', result5);
fprintf('  Exact:  %.10f\n', exact5);
fprintf('  Error:  %.2e\n\n', abs(result5 - exact5));

%% -------------------------------------------------------
% Summary
% -------------------------------------------------------
errors = [abs(result1-exact1), abs(result2-exact2), abs(result3-exact3), ...
          abs(result4-exact4), abs(result5-exact5)];
fprintf('=== Summary ===\n');
fprintf('Max error across all tests: %.2e\n', max(errors));
if max(errors) < 1e-10
    fprintf('All tests PASSED\n');
else
    fprintf('WARNING: some errors exceed 1e-10\n');
end