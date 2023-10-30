function [A, C, S, P, nr] = merge_merged_components(A, C, S, P, nm, neur_id)

%% Function to update fields after merging of components
%
% Originally written by:
% Eftychios A. Pnevmatikakis, Simons Foundation, 2015
%
% ----- Input:
%   A: (p x k matrix, sparse) matrix of spatial footprints (pixel x components)
%   C: (k x n matrix) matrix of temporal footprints (components x frames)
%   S: (k x n matrix) deconvolved activity/spikes (optional)
%   P: (struct) component parameters
%   nm: number of merging operations
%   neur_id: list of updated components
% ----- Output:
%   A: (p x k_new matrix, sparse) updated matrix of spatial footprints (pixel x components)
%   C: (k_new x n matrix) updated matrix of temporal footprints (components x frames)
%   S: (k_new x n matrix) updated deconvolved activity/spikes (optional)
%   P: (struct) updated component parameters
%   nr: new number of components
% ----------

nr = size(A,2);

A = [A(:,1:nr),A_merged,A(:,nr+1:end)];
C = [C(1:nr,:);C_merged;C(nr+1:end,:)];
A(:,neur_id) = [];
C(neur_id,:) = [];

if nargin < 7
    S = [];
    if nargout == 6
        warning('Merged spikes matrix is returned as empty because the original matrix was not provided.');
    end
else
    S = [S(1:nr,:);S_merged];
    S(neur_id,:) = [];
end

if strcmpi(options.deconv_method,'constrained_foopsi') || strcmpi(options.deconv_method,'MCEM_foopsi')
    P.b(neur_id) = [];
    P.b(nr - length(neur_id) + (1:nm)) = P_merged.b;
    P.gn(neur_id) = [];
    P.gn(nr - length(neur_id) + (1:nm)) = P_merged.gn;
    P.c1(neur_id) = [];
    P.c1(nr - length(neur_id) + (1:nm)) = P_merged.c1;
    P.neuron_sn(neur_id) = [];
    P.neuron_sn(nr - length(neur_id) + (1:nm)) = P_merged.neuron_sn;
end
nr = nr - length(neur_id) + nm;

end