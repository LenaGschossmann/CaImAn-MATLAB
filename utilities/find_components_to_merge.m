function [merged_ROIs, nm] = find_components_to_merge(A, C, thr, mx)

%% Function to define components/ROIs that should be merged together
%
% Originally written by:
% Eftychios A. Pnevmatikakis, Simons Foundation, 2015
%
% ----- Input:
%   A: (p x k matrix, sparse) matrix of spatial footprints (pixel x components)
%   C: (k x n matrix) matrix of temporal footprints (components x frames)
%   thr: temporal correlation threshold for merging of components
%   mx: maximum number of merging operations
% ----- Output:
%   merged_ROIs:  list of old components that were merged
%   nm:  number of merging operations
% ----------


nr = size(A,2);

A_corr = triu(A(:,1:nr)'*A(:,1:nr));
A_corr(1:nr+1:nr^2) = 0;
FF2 = A_corr > 0;                                   % find graph of overlapping spatial components

%C_corr = corr(full(C(1:nr,:)'));
C_corr = zeros(nr);
for i = 1:nr
    overlap_indeces = find(A_corr(i,:));
    if ~isempty(overlap_indeces)
        corr_values = corr(C(i,:)',C(overlap_indeces,:)');
        C_corr(i,overlap_indeces) = corr_values;
        C_corr(overlap_indeces,i) = corr_values;
    end
end
FF1 = triu(C_corr)>= thr;                           % find graph of strongly correlated temporal components


FF3 = and(FF1,FF2);                                 % intersect the two graphs
[l,c] = graph_connected_comp(sparse(FF3+FF3'));     % extract connected components
MC = [];
for i = 1:c
    if length(find(l==i))>1
        MC = [MC,(l==i)'];
    end
end

cor = zeros(size(MC,2),1);
for i = 1:length(cor)
    fm = find(MC(:,i));
    for j1 = 1:length(fm)
        for j2 = j1+1:length(fm)
            cor(i) = cor(i) + C_corr(fm(j1),fm(j2));
        end
    end
end

[~,ind] = sort(cor,'descend');
nm = min(length(ind),mx);   % number of merging operations
merged_ROIs = cell(nm,1);
for i = 1:nm
    merged_ROIs{i} = find(MC(:,ind(i)));
end

end