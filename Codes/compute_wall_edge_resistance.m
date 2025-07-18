function [triangle_info, Res] = compute_wall_edge_resistance(nodes_ext, labels_ext, nodes, node_labels, N, eta, y_max)

dt        = delaunayTriangulation(nodes_ext(:,1), nodes_ext(:,2));
triangles = dt.ConnectivityList;
triangle_info = [];

% 3) Filter triangles
for k = 1:size(triangles,1)
    idx = triangles(k,:);
    lbls = labels_ext(idx);
    parts = ceil(abs(lbls)/4);
    is_orig = lbls>0 & lbls<=4*N;
    num_orig = sum(is_orig);
    num_mirr = 3 - num_orig;
    % reject pure-original when N>1
    if N>1 && num_orig==3 && numel(unique(parts))==1
        continue;
    end
    % permit exactly 1 original + 2 same-particle mirrors
    if num_mirr==2 && num_orig==1
        op = parts(is_orig); mp = parts(~is_orig);
        if ~(numel(unique(mp))==1 && mp(1)==op && ...
             (all(lbls(~is_orig)<0) || all(lbls(~is_orig)>4*N)))
            continue;
        end
    elseif num_mirr>0 && num_mirr<3
        % drop triangles with 2 originals? keep only 1orig2mirr or >=2orig
        if num_orig<2, continue; end
    elseif num_mirr==3
        continue;
    end
    % need at least two originals
    if num_orig<2, continue; end
    % check same-particle pair or mirrored match
    parts_o = parts(is_orig);
    accept = false;
    if any(histcounts(parts_o,0.5:1:N+0.5)==2)
        accept = true;
    elseif num_orig==2 && num_mirr==1
        if any(ismember(ceil(lbls(~is_orig)/4), unique(parts_o)))
            accept = true;
        end
    end
    if ~accept, continue; end
    % compute area
    coords = nodes_ext(idx,:);
    A = 0.5 * abs(det([coords(2,:)-coords(1,:); coords(3,:)-coords(1,:)]));
    triangle_info(end+1,:) = [lbls' A];
end

% 4) Initialize resistance
num_nodes = 4*N;
Res = inf(num_nodes, num_nodes);

% 5) Accumulate and assign
for p = 1:N
    n1 = 4*(p-1)+1; n2 = n1+1; n3 = n1+2; n4 = n1+3;
    A_l1=0; A_r1=0; A_l4=0; A_r4=0;
    for t = 1:size(triangle_info,1)
        lbls = triangle_info(t,1:3);
        A    = triangle_info(t,4);
        o_lbls = lbls(lbls>0 & lbls<=4*N);
        parts_o = ceil(o_lbls/4);
        cx = mean(nodes(ismember(node_labels,o_lbls),1));
        % two originals same particle
 % two originals same particle
        if numel(o_lbls)==2 && numel(unique(parts_o))==1
            if all(ismember([n1,n2],o_lbls))
                if nodes(n2,1) < nodes(n1,1), A_l1 = A_l1 + A; else, A_r1 = A_r1 + A; end
            end
            if all(ismember([n1,n3],o_lbls))
                if nodes(n3,1) > nodes(n1,1), A_r1 = A_r1 + A; else, A_l1 = A_l1 + A; end
            end
            if all(ismember([n4,n2],o_lbls))
                if nodes(n2,1) < nodes(n4,1), A_l4 = A_l4 + A; else, A_r4 = A_r4 + A; end
            end
            if all(ismember([n4,n3],o_lbls))
                if nodes(n3,1) > nodes(n4,1), A_r4 = A_r4 + A; else, A_l4 = A_l4 + A; end
            end
        end
        % three originals
        if numel(o_lbls)==3
            if ismember(n1,lbls) && ismember(n2,lbls) && cx<nodes(n1,1), A_l1=A_l1+A; end
            if ismember(n1,lbls) && ismember(n3,lbls) && cx>nodes(n1,1), A_r1=A_r1+A; end
            if ismember(n4,lbls) && ismember(n2,lbls) && cx<nodes(n4,1), A_l4=A_l4+A; end
            if ismember(n4,lbls) && ismember(n3,lbls) && cx>nodes(n4,1), A_r4=A_r4+A; end
        end
    end
    if A_l1>0, R12=8*eta*y_max/(A_l1^2); Res(n1,n2)=R12; Res(n2,n1)=R12; end
    if A_r1>0, R13=8*eta*y_max/(A_r1^2); Res(n1,n3)=R13; Res(n3,n1)=R13; end
    if A_l4>0, R42=8*eta*y_max/(A_l4^2); Res(n4,n2)=R42; Res(n2,n4)=R42; end
    if A_r4>0, R43=8*eta*y_max/(A_r4^2); Res(n4,n3)=R43; Res(n3,n4)=R43; end
end
end
