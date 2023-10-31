function [A_merged, C_merged, P_merged, neur_id] = update_merged_components(Y, A,b, C, f, P, merged_ROIs, nm, options)

%% Function to update spatial and temporal footprints of components/ROIs that are merged together
%
% Originally written by:
% Eftychios A. Pnevmatikakis, Simons Foundation, 2015
%
% ----- Input:
%   Y: (k x n matrix) matrix of raw signal of components (components x frames; only required for fast_merge)
%   A: (p x k matrix, sparse) matrix of spatial footprints (pixel x components)
%   b: (p x 1 matrix) spatial background (pixel x 1)
%   C: (k x n matrix) matrix of temporal footprints (components x frames)
%   f: (1 x n matrix) temporal background (1 x frames)
%   P: (struct) component parameters
%   merged_ROIs: (cell) list of IDs that are merged into new ones
%   nm: number of merging operations
%   options: (struct) algorithm parameters
% ----- Output:
%   A_merged: (p x nm matrix) matrix of merged spatial footprints (pixel x merged components)
%   C_merged: (nm x T matrix) matrix of merged temporal footprints (merged components x frames)
%   S_merged: (nm x T matrix) matrix of merged activity/spikes (merged components x frames)
%   P_merged: (struct) parameters of (new) merged components
%   neur_id: list of updated components
% ----------

d = size(A,1);
T = size(C,2);

A_merged = zeros(d,nm);
C_merged = zeros(nm,T);
S_merged = zeros(nm,T);
Y_merged = zeros(nm,T);
Df_merged = zeros(nm,1);

if strcmpi(options.deconv_method,'constrained_foopsi')
    P_merged.gn = cell(nm,1);
    P_merged.b = cell(nm,1);
    P_merged.c1 = cell(nm,1);
    P_merged.neuron_sn = cell(nm,1);
end

if ~options.fast_merge
    Y_res = Y - A*C;
end


for i = 1:nm
    % merged_ROIs{i} = find(MC(:,ind(i)));
    nC = sqrt(sum(C(merged_ROIs{i},:).^2,2));
    if options.fast_merge
        aa = sum(A(:,merged_ROIs{i})*spdiags(nC,0,length(nC),length(nC)),2);
        for iter = 1:10
            cc = (aa'*A(:,merged_ROIs{i}))*C(merged_ROIs{i},:)/sum(aa.^2);
            aa = A(:,merged_ROIs{i})*(C(merged_ROIs{i},:)*cc')/norm(cc)^2;
        end
        na = sqrt(sum(aa.^2)/max(sum(A(:,merged_ROIs{i}).^2)));
        aa = aa/na;
        %[cc,b_temp,c1_temp,g_temp,sn_temp,ss] = constrained_foopsi(cc);
        cc = na*cc';
        ss = cc;
    else
        A_merged(:,i) = sum(A(:,merged_ROIs{i})*spdiags(nC,0,length(nC),length(nC)),2);
        Y_res = Y_res + A(:,merged_ROIs{i})*C(merged_ROIs{i},:);
        ff = find(A_merged(:,i));
        Pmr = P;
        if isfield(Pmr,'unsaturatedPix')
            px = intersect(Pmr.unsaturatedPix,ff);
            Pmr.unsaturatedPix = zeros(length(px),1);
            for pxi = 1:length(px)
                Pmr.unsaturatedPix(pxi) = find(ff == px(pxi));
            end
        end
        cc = update_temporal_components(Y_res(ff,:),A_merged(ff,i),b(ff,:),median(spdiags(nC,0,length(nC),length(nC))\C(merged_ROIs{i},:)),f,Pmr,options);
        [aa,bb] = update_spatial_components(Y_res,cc,f,A_merged(:,i),P,options);
        [cc,~,Ptemp,ss] = update_temporal_components(Y_res(ff,:),aa(ff),bb(ff,:),cc,f,Pmr,options);
    end
    A_merged(:,i) = aa;
    C_merged(i,:) = cc;
    S_merged(i,:) = ss;
    if options.fast_merge
        Y_merged(i,:) = yy;
        Df_merged(i) = df;
    end

    if strcmpi(options.deconv_method,'constrained_foopsi') || strcmpi(options.deconv_method,'MCEM_foopsi')
        if options.fast_merge
            P_merged.gn{i} = 0; %g_temp;   % do not perform deconvolution during merging
            P_merged.b{i} = 0;  %b_temp;
            P_merged.c1{i} = 0; %c1_temp;
            P_merged.neuron_sn{i} = 0; %sn_temp;
        else
            P_merged.gn{i} = Ptemp.gn{1};
            P_merged.b{i} = Ptemp.b{1};
            P_merged.c1{i} = Ptemp.c1{1};
            P_merged.neuron_sn{i} = Ptemp.neuron_sn{1};
            if i < nm
                Y_res(ff,:) = Y_res(ff,:) - aa(ff)*cc;
            end
        end
    end
end


neur_id = unique(cell2mat(merged_ROIs));


end