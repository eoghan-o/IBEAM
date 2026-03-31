function [g, un] = build_graph_and_normals(vertices, local_face)
    n_v     = size(vertices, 1);
    n_faces = size(local_face, 1);
    poly_centroid = mean(vertices, 1);

    % Compute outward unit normals
    un = zeros(n_faces, 3);
    for f = 1:n_faces
        idx    = local_face(f,:);
        vi     = vertices(idx(1),:);
        vj     = vertices(idx(2),:);
        vk     = vertices(idx(3),:);
        normal = cross(vj - vi, vk - vi);
        fc     = mean(vertices(idx,:), 1);
        if dot(normal, fc - poly_centroid) < 0
            normal = -normal;
        end
        un(f,:) = normal / norm(normal);
    end

    % Build g{i}: incident faces in CCW order around vertex i as seen from outside
    % Strategy: chain faces by shared edges, then check orientation
    g = cell(n_v, 1);
    for i = 1:n_v
        % Find incident faces
        incident = [];
        for f = 1:n_faces
            if ismember(i, local_face(f,:))
                incident(end+1) = f;
            end
        end
        n_inc = length(incident);

        if n_inc <= 1
            g{i} = incident;
            continue
        end

        % Build edge-sharing adjacency between incident faces around vertex i
        % Two faces are adjacent around i if they share exactly one other vertex
        neighbors = zeros(n_inc, 2);
        nb_count  = zeros(n_inc, 1);
        shared_verts = zeros(n_inc, n_inc);  % shared_verts(a,b) = shared vertex between face a and b
        for fi = 1:n_inc
            for fj = fi+1:n_inc
                fa     = local_face(incident(fi),:);
                fb     = local_face(incident(fj),:);
                shared = setdiff(intersect(fa, fb), i);
                if length(shared) == 1
                    nb_count(fi) = nb_count(fi) + 1;
                    nb_count(fj) = nb_count(fj) + 1;
                    neighbors(fi, nb_count(fi)) = fj;
                    neighbors(fj, nb_count(fj)) = fi;
                    shared_verts(fi,fj) = shared;
                    shared_verts(fj,fi) = shared;
                end
            end
        end

        % Walk the ring starting from face 1
        order      = zeros(1, n_inc);
        order(1)   = 1;
        visited    = false(n_inc, 1);
        visited(1) = true;
        prev       = 0;  % no previous face at start
        for step = 2:n_inc
            curr = order(step-1);
            % Pick the neighbor that isn't where we came from
            next = 0;
            for nb = 1:2
                cand = neighbors(curr, nb);
                if cand > 0 && ~visited(cand)
                    next = cand;
                    break;
                end
            end
            if next == 0
                remaining  = find(~visited);
                next       = remaining(1);
            end
            order(step) = next;
            visited(next) = true;
        end
        sorted_faces = incident(order);

        % Check CCW orientation of the ring around vertex i from outside:
        % Take first two faces in ring. Their shared edge goes from vertex i
        % to shared_v. The cross product (shared_v - i) x (next_v - i)
        % should point outward from the polyhedron.
        f1       = order(1);
        f2       = order(2);
        sv       = shared_verts(f1, f2);  % shared vertex between face 1 and 2
        vi       = vertices(i,:);
        out_dir  = vi - poly_centroid;

        if sv > 0
            % Get a vertex in face 1 that is NOT i and NOT sv
            fa      = local_face(incident(f1),:);
            other_v = fa(fa ~= i & fa ~= sv);
            if ~isempty(other_v)
                e1    = vertices(sv,        :) - vi;
                e2    = vertices(other_v(1),:) - vi;
                cross_e = cross(e1, e2);
                if dot(cross_e, out_dir) < 0
                    % Wrong orientation — reverse the ring
                    sorted_faces = fliplr(sorted_faces);
                end
            end
        end

        g{i} = sorted_faces;
    end
end