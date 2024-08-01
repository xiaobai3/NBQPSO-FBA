function [index, model_bio] = preprocessModel(model, biomass_rxn, target_rxn, bound_glu, bound_O2, bound_ATPM, glc_rxn, o2_rxn)
% Generates an index vector for EBQPSO algorithm by removing essential reactions and genes
% Input:
%   - model: COBRA stoichiometric model in standard structure form
%   - biomass_rxn: Biomass reaction
%   - target_rxn: Target reaction
%   - bound_glu: Input flux value for glucose
%   - bound_O2: Input flux value for oxygen
%   - bound_ATPM: Input flux value for maintenance ATP
% Output:
%   - index: Index vector for EBQPSO algorithm
    
    warning off optim:linprog:WillRunDiffAlg
    
    model_bio = setBounds(model, glc_rxn, o2_rxn, bound_glu, bound_O2, bound_ATPM);
    essential_reactions = getEssentialReactions(model_bio);
    
    % Remove essential, ATPM, target, and biomass reactions
    rxns_to_remove = union({biomass_rxn, 'ATPM', target_rxn}, model_bio.rxns(essential_reactions==1));
    [~, rxns_to_remove_ids] = ismember(rxns_to_remove, model_bio.rxns);
    ind = setdiff(1:length(model_bio.rxns), rxns_to_remove_ids);
    
    % Remove reactions not associated with genes and extracellular exchange reactions
    ind = removeNoGeneReactions(model_bio, ind);
    ind = removeExtracellularExchangeReactions(model_bio, ind);
    
    index = ind';
end

function model = setBounds(model, glc_rxn, o2_rxn, bound_glu, bound_O2, bound_ATPM)
% Sets bounds for glucose, oxygen, and ATPM reactions
% Input:
%   - model: COBRA stoichiometric model
%   - bound_glu: Input flux value for glucose
%   - bound_O2: Input flux value for oxygen
%   - bound_ATPM: Input flux value for maintenance ATP
% Output:
%   - model: COBRA model with updated bounds

    model = changeRxnBounds(model, glc_rxn, bound_glu, 'l');
    model = changeRxnBounds(model, o2_rxn, bound_O2, 'l');
    model = changeRxnBounds(model, 'ATPM', bound_ATPM, 'l');
end

function essential_reactions = getEssentialReactions(model)
% Identifies essential reactions based on FBA solution
% Input:
%   - model: COBRA stoichiometric model
% Output:
%   - essential_reactions: Binary vector indicating essential reactions

    essential_reactions = zeros(length(model.rxns), 1);
    for i = 1:length(model.rxns)
        model_del = model;
        model_del = changeRxnBounds(model_del, model.rxns{i}, 0, 'b');
        FBA_solution = optimizeCbModel(model_del, 'max');
        if FBA_solution.f == 0 || isnan(FBA_solution.f)
            essential_reactions(i) = 1; % Essential reaction
        end
    end
end


function ind = removeNoGeneReactions(model, ind)
% Removes reactions not associated with genes from the index list
% Input:
%   - model: COBRA stoichiometric model
%   - ind: Index list of reactions
% Output:
%   - ind: Updated index list with reactions not associated with genes removed

    no_gene_indices = find(cellfun(@isempty, model.grRules));
    ind = setdiff(ind, no_gene_indices);
end

function ind = removeExtracellularExchangeReactions(model, ind)
% Removes extracellular exchange reactions from the index list
% Input:
%   - model: COBRA stoichiometric model
%   - ind: Index list of reactions
% Output:
%   - ind: Updated index list with extracellular exchange reactions removed

    extracellular_rxn_indices = find(strcmp(model.subSystems, 'Extracellular exchange'));
    ind = setdiff(ind, extracellular_rxn_indices);
end

