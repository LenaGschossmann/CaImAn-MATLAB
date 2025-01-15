function [A,C,nr, merged_ROIs, P,S, neur_id] = merge_components(Y,A,b,C,f,P,S,options,merged_ROIs)

% merging of spatially overlapping components that have highly correlated tmeporal activity
% The correlation threshold for merging overlapping components is user specified in P.merge_thr (default value 0.85)
% Inputs:
% Y:            raw data
% A:            matrix of spatial components
% b:            spatial background
% C:            matrix of temporal components
% f:            temporal background
% P:            struct for neuron parameters
% S:            deconvolved activity/spikes (optional)
% options:      struct for algorithm parameters
% merged_ROIs:  list of IDs that are merged into new ones

% Outputs:
% A:            matrix of new spatial components
% C:            matrix of new temporal components
% nr:           new number of components
% merged_ROIs:  list of IDs that are merged into new ones
% P:            new parameter struct
% S:            matrix of new deconvolved/activity spikes
% neur_id:      list of updated components
% Adapted from code written by:
% Eftychios A. Pnevmatikakis, Simons Foundation, 2015

defoptions = CNMFSetParms;
if isempty(options), options = []; end
if ~isfield(options,'d1') || isempty(options.d1); d1 = input('What is the total number of rows? \n'); else d1 = options.d1; end          % # of rows
if ~isfield(options,'d2') || isempty(options.d2); d2 = input('What is the total number of columns? \n'); else d2 = options.d2; end       % # of columns
if ~isfield(options,'merge_thr') || isempty(options.merge_thr); thr = defoptions.merge_thr; else thr = options.merge_thr; end     % merging threshold
if ~isfield(options,'max_merg'); mx = 50; else mx = options.max_merg; end           % maximum merging operations
if ~isfield(options,'deconv_method') || isempty(options.deconv_method); options.deconv_method = defoptions.deconv_method; end
if ~isfield(options,'fast_merge') || isempty(options.fast_merge); options.fast_merge = defoptions.fast_merge; end  % flag for using fast merging

nr = size(A,2);

if nr == 0
    merged_ROIs = [];
    return
end

if isempty(merged_ROIs)
    [merged_ROIs, nm] = find_components_to_merge();
else % merged_ROIs is provided, allowing for custom defining merged_ROIs.
    nm = length(merged_ROIs);
end

% Update ROIs that are merged together
[A_merged, C_merged, S_merged, P_merged, neur_id] = update_merged_components(Y, A, b, C, f, P, merged_ROIs, nm, options);

% Merge the newly-merged components back with the not merged ones
[A, C, S, P] = add_merged_components(A, A_merged, C, C_merged, S, S_merged, P, P_merged, nm, neur_id, options);

nr = size(C,1);

