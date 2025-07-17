function edges_ascending = getEdgesAscending(edges,nodes)
% Prepare edge list such that i is lower y than j
edges_ascending = zeros(size(edges));
for k = 1:size(edges, 1)
    i = edges(k, 1);
    j = edges(k, 2);
    
    if nodes(i, 2) <= nodes(j, 2)
        edges_ascending(k, :) = [i, j];
    else
        edges_ascending(k, :) = [j, i];
    end
end

% Now sort rows based on y-value of node 'i'
y_i = nodes(edges_ascending(:,1), 2); % y-coordinates of all i's
[~, sort_idx] = sort(y_i);            % indices that would sort y_i ascending
edges_ascending = edges_ascending(sort_idx, :); % apply the sort