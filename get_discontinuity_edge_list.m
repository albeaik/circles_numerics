function output_edge_list = get_discontinuity_edge_list(DT, density, maxdtdensity)
% DT generated from matlab built in functions
% density is a list of density values per triangle
% discontinuity assumption: delta(density) > maxdtdensity ==> discont

    %discontinuity assumption: delta(density) > maxdtdensity ==> discont
    %maxdtdensity = 0.5;

    %identify searchable trig adjacancy list
    trig_neighbor_sets = [[1:size(DT, 1)]', neighbors(DT)]; %identify neighboring trigs
    trig_neighbor_list = [trig_neighbor_sets(:, [1, 2]); trig_neighbor_sets(:, [1, 3]); trig_neighbor_sets(:, [1, 4])]; %turn into trig-trig ID list
    trig_neighbor_list_clean = trig_neighbor_list(sum(isnan(trig_neighbor_list),2) == 0, :); %remove boundary elements
    trig_neighbor_list_sorted = sort(trig_neighbor_list_clean, 2);  %sort to simplify search

    unique_trig_neighbor_list = unique(trig_neighbor_list_sorted, 'rows');  %idenitfy unique sets of neighbors

    %identify discontinuities
    discontinuity_flags = zeros([size(unique_trig_neighbor_list, 1), 1]);
    for i = 1:size(unique_trig_neighbor_list, 1)
        density_trig_1 = density(unique_trig_neighbor_list(i, 1));
        density_trig_2 = density(unique_trig_neighbor_list(i, 2));

        if(abs(density_trig_1 - density_trig_2) > maxdtdensity)
            discontinuity_flags(i) = 1;
        end
    end

    %extract edges of discontinuities
    unique_trig_neighbor_discont_list = unique_trig_neighbor_list(discontinuity_flags == 1, :);
    discontinuity_edge_list = zeros(size(unique_trig_neighbor_discont_list));
    for i = 1:size(unique_trig_neighbor_discont_list, 1)
        trig_1 = DT(unique_trig_neighbor_discont_list(i, 1), :);
        trig_2 = DT(unique_trig_neighbor_discont_list(i, 2), :);
        discontinuity_edge_list(i, :) = intersect(trig_1, trig_2);
    end

    output_edge_list = discontinuity_edge_list;
    
end
