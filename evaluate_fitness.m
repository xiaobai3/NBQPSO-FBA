function [Fitvalue, bioValue, ProdValue] = evaluate_fitness(indexSelected, bound_biom)
% FITNESSFUN_EBQPSO_RXNS Calculates the fitness value for EBQPSO algorithm
% Input:
%   - indexSelected: Binary matrix representing selected reactions for deletion
% Output:
%   - Fitvalue: Fitness value for each individual
%   - bioValue: Biomass production for each individual
%   - ProdValue: Target product formation for each individual

    % Access global variables
    global model_bio;
    global bound_glu;
    global target_rxn;
    global toler;
    global indexing;
    
    % Initialize variables
    popsize = size(indexSelected, 1);
    bioValue = zeros(popsize, 1);
    ProdValue = zeros(popsize, 1);
    Fitvalue = zeros(popsize, 1);
    
    % Loop through each individual in the population
    for i = 1:popsize
        % Get reaction deletions for the current individual
        deletions = model_bio.rxns(indexing(indexSelected(i,:)));
        nDel = length(deletions);
        model_KO = model_bio;

        % Apply reaction deletions to the model
        for j = 1:nDel
            model_KO = changeRxnBounds(model_KO, deletions{j}, 0, 'b');
        end

        % Optimize the biomass production in the knockout model
        AfterKO = optimizeCbModel(model_KO);
        growthRate = roundn(AfterKO.f, -4);
        if isnan(growthRate)
            growthRate = 0;
        end
        bioValue(i) = growthRate;
        
        % If the growth rate is above a threshold, optimize the target product formation
        if (growthRate > bound_biom)
            round_off = floor(AfterKO.f / toler) * toler;
            model_KO = changeRxnBounds(model_KO, model_KO.rxns(model_KO.c == 1), round_off, 'l');
            model_KO = changeObjective(model_KO, target_rxn);
            Max_Sol = optimizeCbModel(model_KO, 'max');
            Prod_Max = Max_Sol.f;
        else
            Prod_Max = 0;
        end
        
        % Calculate and store the production value and fitness value for each individual
        ProdValue(i) = roundn(Prod_Max, -4);
%         Fitvalue(i) = roundn(bioValue(i) * ProdValue(i) / abs(bound_glu), -4);
        Fitvalue(i) = ProdValue(i);
    end
end