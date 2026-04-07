function export_to_paraview(nodes, elements, u_h, u_ex, filename)
% EXPORT_TO_PARAVIEW  Write FEM solution to VTK Unstructured Grid (.vtu)
%   for visualization in ParaView.
%
%   nodes    - (n_nodes x 3) array of node coordinates [x, y, z]
%   elements - (n_elems x 4) array of tet element connectivity (1-indexed)
%   u_h      - (n_nodes x 1) FEM solution vector
%   u_ex     - (n_nodes x 1) exact solution at nodes
%   filename - output filename, e.g. 'solution.vtu'
%
%   Point data written:  u_h, u_exact, pointwise_error
%
%   Usage:
%     export_to_paraview(nodes, elements, u_h, u_ex, 'solution_k5.vtu');

    if nargin < 5, filename = 'solution.vtu'; end
    if ~endsWith(filename, '.vtu'), filename = [filename '.vtu']; end

    n_nodes = size(nodes, 1);
    n_elems = size(elements, 1);
    u_err   = abs(u_h - u_ex);

    % VTK cell type 10 = VTK_TETRA (4-node tetrahedron)
    VTK_TETRA = 10;

    fid = fopen(filename, 'w');
    if fid == -1
        error('Cannot open file: %s', filename);
    end

    %% Header
    fprintf(fid, '<?xml version="1.0"?>\n');
    fprintf(fid, '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">\n');
    fprintf(fid, '  <UnstructuredGrid>\n');
    fprintf(fid, '    <Piece NumberOfPoints="%d" NumberOfCells="%d">\n', n_nodes, n_elems);

    %% Points
    fprintf(fid, '      <Points>\n');
    fprintf(fid, '        <DataArray type="Float64" NumberOfComponents="3" format="ascii">\n');
    fprintf(fid, '          ');
    fprintf(fid, '%f %f %f\n', nodes');
    fprintf(fid, '        </DataArray>\n');
    fprintf(fid, '      </Points>\n');

    %% Cells (connectivity, offsets, types)
    % VTK uses 0-based indexing
    conn = (elements - 1)';  % 4 x n_elems, 0-indexed

    fprintf(fid, '      <Cells>\n');

    % Connectivity
    fprintf(fid, '        <DataArray type="Int32" Name="connectivity" format="ascii">\n');
    fprintf(fid, '          ');
    fprintf(fid, '%d %d %d %d\n', conn);
    fprintf(fid, '        </DataArray>\n');

    % Offsets (cumulative node count per cell)
    fprintf(fid, '        <DataArray type="Int32" Name="offsets" format="ascii">\n');
    fprintf(fid, '          ');
    fprintf(fid, '%d\n', (4:4:4*n_elems)');
    fprintf(fid, '        </DataArray>\n');

    % Cell types
    fprintf(fid, '        <DataArray type="UInt8" Name="types" format="ascii">\n');
    fprintf(fid, '          ');
    fprintf(fid, '%d\n', repmat(VTK_TETRA, n_elems, 1));
    fprintf(fid, '        </DataArray>\n');

    fprintf(fid, '      </Cells>\n');

    %% Point data (nodal scalar fields)
    fprintf(fid, '      <PointData Scalars="u_h">\n');

    write_scalar(fid, 'u_h',            u_h,    n_nodes);
    write_scalar(fid, 'u_exact',        u_ex,   n_nodes);
    write_scalar(fid, 'pointwise_error', u_err, n_nodes);

    fprintf(fid, '      </PointData>\n');

    %% Footer
    fprintf(fid, '    </Piece>\n');
    fprintf(fid, '  </UnstructuredGrid>\n');
    fprintf(fid, '</VTKFile>\n');

    fclose(fid);
    fprintf('Exported to %s  (%d nodes, %d tets)\n', filename, n_nodes, n_elems);
end

%% Helper: write a named nodal scalar array
function write_scalar(fid, name, data, n_nodes)
    fprintf(fid, '        <DataArray type="Float64" Name="%s" format="ascii" NumberOfComponents="1">\n', name);
    fprintf(fid, '          ');
    fprintf(fid, '%.10e\n', data);
    fprintf(fid, '        </DataArray>\n');
end

function tf = endsWith(str, suffix)
    tf = length(str) >= length(suffix) && strcmp(str(end-length(suffix)+1:end), suffix);
end