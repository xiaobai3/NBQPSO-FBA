clc; clear; close all;

warning off optim:linprog:WillRunDiffAlg
%% Set the parameters
% Initialize the COBRA TOOLbox if not initialized
if ~exist('initCobraToolbox', 'file')
    initCobraToolbox(false);
end
changeCobraSolver('gurobi');

global model_bio;
global bound_glu;
global biomass_rxn;
global target_rxn;
global toler;
global indexing;
global ind_bio;
global ind_tar;

% Load the COBRA model and set parameters
% load('iJR904.mat');
% model = iJR904;
% model=readCbModel('iJR904.xml');
load('iAF1260.mat');
model = iAF1260;

% load('iML1515.mat')

% % iJR904 model
% biomass_rxn = 'BIOMASS_Ecoli'; % biomass in iJR904 model
% target_rxn = 'EX_lac__D_e'; % D-lactate as a target reaction in iJR904
% model

% iAF1260 model
biomass_rxn = 'BIOMASS_Ec_iAF1260_core_59p81M'; % biomass in iAF1260 model
target_rxn = 'EX_ac_e'; % iAF1260中以acetate为目标反应
% target_rxn = 'EX_succ_e'; % iAF1260中以Succinate为目标反应

% % iML1515 model
% % target_rxn='EX_succ_e';
% target_rxn='EX_ac_e';
% biomass_rxn=model.rxns{model.c==1};

[~, ind_bio] = ismember(biomass_rxn, model.rxns); % ind_bio is the rxn number of biomass_rxn
[~, ind_tar] = ismember(target_rxn, model.rxns); % ind_tar is the rxn number of target_rxn

% iAF1260 MODEL
glc_rxn = 'EX_glc__D_e'; % iAF1260
o2_rxn = 'EX_o2_e'; % iAF1260

% % iML1515 model
% o2_rxn='EX_o2_e';
% glc_rxn='EX_glc__D_e';

bound_glu = -10;
bound_O2 = -18.5; 
% bound_O2 = -100; 
bound_biom = 0.1;
% bound_ATPM = 7.6; % iJR904 ATPM
bound_ATPM = 8.39; %  iAF1260 ATPM
% bound_ATPM = 6.86; % iML1515 ATPM
[indexing, model_bio] = preprocessModel(model, biomass_rxn, target_rxn, bound_glu, bound_O2, bound_ATPM, bound_biom, glc_rxn, o2_rxn);
% [~, model_bio] = preprocessModel(model, biomass_rxn, target_rxn, bound_glu, bound_O2, bound_ATPM, glc_rxn, o2_rxn);
% [~, model_bio] = preprocessModel(model, biomass_rxn, target_rxn, bound_glu, bound_O2, bound_ATPM, bound_biom, glc_rxn, o2_rxn)
% indexing=xlsread('iAF1260_succinate_iter6');
% load('indexing_iML1515_succinate_myPreProM.mat'); % load indexing
% load('iML1515_preProssedModel.mat');% load model_bio
% load('iML1515_preProssedModel_100.mat');% load model_bio
% load('iML1515_acetate_PreProM_100.mat');
% load('indexing_iML1515_acetate_PreProM.mat');

% model_bio = changeRxnBounds(model, glc_rxn, bound_glu, 'b');
% model_bio = changeRxnBounds(model_bio, biomass_rxn, bound_biom, 'l'); 

% % limit reaction rate in realistic range
% model_bio.lb(model_bio.lb<-100)=-100;
% model_bio.ub(model_bio.ub>100)=100;

% 
% load('candidate_OptDesign.mat');
% candidateRxns = candidate.rxns;
% 
%  % Initialize an empty cell array to store the split strings
% SepCandidateRxns = {};
% 
% % Loop through the input cell array and split strings by '/'
% for i = 1:numel(candidateRxns)
%     % Check if the current element is a string
%     if ischar(candidateRxns{i})
%         % Split the current string by '/'
%         splitStrings = strsplit(candidateRxns{i}, '/');
%         
%         % Add the split strings to the output cell array
%         SepCandidateRxns = [SepCandidateRxns, splitStrings];
%     end
% end
% 
% [~, indexing] = ismember(SepCandidateRxns, model_bio.rxns);

% Parameters
max_fitness_num = 30000;
% popSize = 100;
popSize_list = (100);
% popSize_list = [100,200];
% Itermax = 300;
% popSize = 500;
% beta = 1.4; % The optimal alpha for others
% beta = 1; % The optimal alpha for iAF1260__EX_succ_e
% beta_list = [1.2, 1.4, 1.6]; 
beta_list = (1.4); 
toler = 1e-7;
neighborhoodSize = 5;
% neighborhoodSize_list = [3,5,7,9];

rxnNums = length(indexing);
rxnBit = length(dec2bin(rxnNums));
maxDeletedNum = 5;
dim = rxnBit * maxDeletedNum;
run_times = 20;
fitness_eval = 'evaluate_fitness';

tic
% Main loop
for popSize = popSize_list
    fprintf('popSize: %d\n',popSize);
    Itermax = max_fitness_num/popSize;
    rep.best_fits = zeros(Itermax, run_times);
    rep.best_bios = zeros(Itermax, run_times);
    rep.best_prods = zeros(Itermax, run_times);
    rep.gbest_fits = zeros(run_times, 1);
    rep.gbest_inds = zeros(run_times, maxDeletedNum);
    rep.gbest_bios = zeros(run_times, 1);
    rep.gbest_prods = zeros(run_times, 1);

    % gbest_fits_list = (9.2781);
    % average_fits_list = (9.2545);

    gbest_fits_list = [];
    average_fits_list = [];
    for beta = beta_list
        fprintf('alpha: %d\n',beta);
        for run_t = 1:run_times
            % Initialize particles
            [particles, indexSelected] = Initialize(popSize, rxnBit, maxDeletedNum, rxnNums);

            % Calculate fitness
            [fitness, bioValue, ProdValue] = evaluate_fitness(indexSelected, bound_biom);

            % Initialize pbest and gbest
            [pbest, pbest_ind, pbest_fit, pbest_bio, pbest_prod, gbest, gbest_ind, gbest_fit, gbest_bio, gbest_prod] = initialize_pbest_gbest(particles, indexSelected, fitness, bioValue, ProdValue);

            best_fitness=zeros(Itermax,1);
            best_bio = zeros(Itermax,1);
            best_prod = zeros(Itermax,1);
            best_ind = zeros(Itermax,maxDeletedNum);

            for t = 1:Itermax
                % Calculate the neighbor best for each particle
                nbest = CalculateNBest(particles, popSize, pbest_fit, neighborhoodSize);

                % Calculate the mbest
                mbest = Get_mbest(pbest);

                % Update particle positions
                [particles, indexSelected] = UpdateParticles(popSize, maxDeletedNum, rxnBit, rxnNums, particles, indexSelected, pbest, pbest_fit, gbest, nbest, mbest,beta);

                % Calculate fitness for updated particles
                [fitness, bioValue, ProdValue] = evaluate_fitness(indexSelected, bound_biom);

                % Update pbest and gbest
                [pbest, pbest_ind, pbest_fit, pbest_bio, pbest_prod, gbest, gbest_ind, gbest_fit, gbest_bio, gbest_prod] = update_pbest_gbest(particles, indexSelected, fitness, bioValue, ProdValue, pbest, pbest_ind, pbest_fit, pbest_bio, pbest_prod, gbest, gbest_ind, gbest_fit, gbest_bio, gbest_prod);

                % Save the results
                best_fitness(t,1) = gbest_fit;    
                best_bio(t) = gbest_bio;
                best_prod(t) = gbest_prod;
                best_ind(t,:) = gbest_ind;

                % Display the results
                if rem(t,50) == 0
                    fprintf('Iter: %f\n', t);
                fprintf('Best fitness: %f\n', gbest_fit);
                fprintf('Best bio: %f\n', gbest_bio);
                fprintf('Best prod: %f\n', gbest_prod);
                end
            end

            % Save the results
            rep.best_fits(:,run_t) = best_fitness;
            rep.best_bios(:,run_t) = best_bio;
            rep.best_prods(:,run_t) = best_prod;

            rep.gbest_fits(run_t) = gbest_fit;
            rep.gbest_inds(run_t,:) = gbest_ind;   
            rep.gbest_bios(run_t,:) = gbest_bio;
            rep.gbest_prods(run_t,:) = gbest_prod;

            disp(['Run：',num2str(run_t)])
            BestdeletedGenes = model_bio.rxns(indexing(gbest_ind));
            rep.gbest_knockouts(run_t,:) = BestdeletedGenes;
            
            fprintf('Optimal bpcy: %f\n',gbest_fit);
            disp('Corresponding reaction deletions:');
            disp(BestdeletedGenes);
        end
        % Figures
        figure()
        % plot(1:Itermax,best_fitness,'-');
        rep.mean_best_fits = mean(rep.best_fits,2);
        plot(1:Itermax,rep.mean_best_fits,'-');
        % Adding Axis Names
        xlabel('Iteration') 
        ylabel('Yield Rate')
        % 添加箭头到坐标轴
        annotation('arrow', [0.13, 0.95], [0.11, 0.11]) % x-Axis arrow
        annotation('arrow', [0.13, 0.13], [0.11, 0.99]) % y-axis arrow
        box off;
    %     output_filename = sprintf('%s__%s_%s%s_%s%s_%s%s',model.description,target_rxn, num2str(Itermax), 'iters',num2str(run_times),'runs', num2str(popSize),'pop');
    %     output_filename = sprintf('%s__%s_%s%s_%s%s_%s%s',model.description,target_rxn, num2str(Itermax), 'iters',num2str(beta * 10),'alpha', num2str(popSize),'pop');
         output_filename = sprintf('BQPSO_%s__%s_%s%s_%s%s_%s%s_%s%s__original',model.description,target_rxn, num2str(popSize),'pop', num2str(Itermax), 'iters',num2str(beta * 10),'alpha', num2str(neighborhoodSize),'neighborhoodSize');
        % saveas(gcf, 'iAF1260_acetate_300iters_10runs', 'png');
        % saveas(gcf, 'iAF1260_acetate_300iters_10runs', 'fig');
        saveas(gcf, output_filename, 'emf');
        saveas(gcf, output_filename, 'fig');
        grid on;
        grid minor;

        output_filename_mat = sprintf('%s%s',output_filename,'.mat');
        save(output_filename_mat, 'rep');

        [ggbest_fit,ind]= max(rep.gbest_fits);
        ggbest_ind = rep.gbest_inds(ind,:);
        ggbest_bio = rep.gbest_bios(ind);
        ggbest_prod = rep.gbest_prods(ind);
        Bestdeletions = model_bio.rxns(indexing(ggbest_ind));
        fprintf('popSize: %d\n', popSize);
        fprintf('Optimal BPCY: %f\n',ggbest_fit);
        fprintf('Corresponding production: %f\n',ggbest_prod);
        fprintf('Corresponding biomass: %f\n',ggbest_bio);
        disp('Corresponding reaction deletions:');
        disp(Bestdeletions);
        mean_gbest_fits = mean(rep.gbest_fits);
        fprintf('Average optimal BPCY: %f\n',mean_gbest_fits);
        fprintf('Average optimal production: %f\n',mean(rep.gbest_prods));
        fprintf('Average optimal biomass: %f\n',mean(rep.gbest_bios));

        gbest_fits_list = [gbest_fits_list,ggbest_fit];
        average_fits_list = [average_fits_list, mean_gbest_fits];
    end

end
executeTime = toc;
disp( ['Average runtime per run: ',num2str(executeTime / run_times)] );
disp( ['Total runtime: ',num2str(executeTime)] );


function [particles, indexSelected] = Initialize(particleNum, rxnBit, maxDeletedNum, rxnNums)
% The Initialize function is responsible for creating the initial particle
% positions (in binary format) and the corresponding indexSelected matrix
% (in decimal format). It takes four input arguments:
% 
% particleNum: The number of particles in the swarm.
% rxnBit: The number of bits required to represent a reaction index in binary format.
% maxDeletedNum: The maximum number of reactions that can be deleted.
% rxnNums: The total number of reactions in the model.

% The function initializes two matrices:
% 
% particles: A matrix with dimensions particleNum x (rxnBit *
% maxDeletedNum), where each row represents a particle and each group of
% rxnBit columns represents a reaction index in binary format.

% indexSelected: A matrix with dimensions particleNum x maxDeletedNum,
% where each row represents a particle and each column represents a
% reaction index in decimal format.

% The function iterates through each particle and reaction index,
% generating a random reaction index (in decimal format) and converting it
% to binary format. The binary representation is then stored in the
% particles matrix, while the decimal representation is stored in the
% indexSelected matrix.

    indexSelected = zeros(particleNum, maxDeletedNum);
    particles = zeros(particleNum, maxDeletedNum * rxnBit);

    for i = 1:particleNum
        temp = randperm(rxnNums, maxDeletedNum);
        indexSelected(i, :) = temp;

        for j = 1:maxDeletedNum
            tmp = dec2bin(indexSelected(i, j), rxnBit) - '0';
            particles(i, (j - 1) * rxnBit + (1:rxnBit)) = tmp;
        end
    end
end
% This code now generates the binary representation for each value in the
% indexSelected matrix and assigns them to the corresponding row in the
% particles matrix in the correct format.

function [pbest, pbest_ind, pbest_fit, pbest_bio, pbest_prod, gbest, gbest_ind, gbest_fit, gbest_bio, gbest_prod] = initialize_pbest_gbest(particles, indexSelected, fitness, bioValue, ProdValue)
    % Initialization of pbest and gbest
    pbest = particles;
    pbest_ind = indexSelected;
    pbest_fit = fitness;
    pbest_bio = bioValue;
    pbest_prod = ProdValue;
    [gbest_fit, ind] = max(fitness);
    gbest = particles(ind, :);
    gbest_ind = indexSelected(ind, :);
    gbest_bio = bioValue(ind);
    gbest_prod = ProdValue(ind);
end

function [particles, indexSelected] = UpdateParticles(particleNum, maxDeletedNum, rxnBit, rxnNums, particles, indexSelected, pbest, pbest_fit, gbest, nbest, mbest, beta)
    % UpdateParticles function updates particle positions and corresponding
    % indexSelected matrix using the provided pbest, nbest, and mbest.
    %
    % particleNum: The number of particles in the swarm.
    % maxDeletedNum: The maximum number of reactions that can be deleted.
    % rxnBit: The number of bits required to represent a reaction index in binary format.
    % rxnNums: The total number of reactions in the model.
    % particles: The current particle positions.
    % indexSelected: The current indexSelected matrix.
    % pbest: The personal best positions for each particle.
    % pbest_fit: The fitness values for pbest.
    % nbest: The neighbor best positions for each particle.
    % mbest: The mean best position.
    % beta: The beta parameter used for calculating the probability of updating a position.
    %
    % The function returns the updated particle positions and corresponding indexSelected matrix.
    P = zeros(particleNum, maxDeletedNum * rxnBit);
    for i = 1:particleNum
        k = randi(particleNum); % Randomly select an integer as another particle for Get_P
        while k == i
            k = randi(particleNum); % Randomly select a non-i integer as another particle for Get_P
        end
        if pbest_fit(k) > pbest_fit(i) % If the selected k-th particle's pbest is better than pbest i
            G = pbest(k,:);
%             G = nbest(k,:);
        else    
            G = gbest;
        end

        % Perform the Get_P2 operation on pbest, G, and nbest to obtain P
%         P(i,:) = Get_P2(pbest(i,:), G, nbest(i,:));
        P(i,:) = Get_P(pbest(i,:), gbest);  % Obtain P by performing crossover on pbest and gbest (for BQPSO case)

        % Update particle positions and corresponding indexSelected values
        for j = 1:maxDeletedNum
            id_sta = (j - 1) * rxnBit + 1;
            id_end = j * rxnBit;

            c = [particles(i,id_sta:id_end); mbest(id_sta:id_end)];
            b = beta * pdist(c, 'hamming') * rxnBit * log(1/rand); 
            Pr = b/rxnBit;
            if Pr > 1
                Pr = 1;
            end

            particles(i,id_sta:id_end) = Transf(P(i,id_sta:id_end),Pr);
            indexSelected(i,j) = bin2dec(num2str(particles(i,id_sta:id_end)));
            if indexSelected(i,j) > rxnNums || indexSelected(i,j) == 0
                indexSelected(i,j) = randperm(rxnNums, 1); % If the updated position value is greater than rxnNums or zero, randomly assign a value from 1 to rxnNums
                tmp = dec2bin(indexSelected(i,j), rxnBit);
                k = id_sta;
                for kk = 1:rxnBit
                    particles(i, k) = str2num(tmp(kk));
                    k = k + 1;
                end
            end
        end
    end
end


function [pbest, pbest_ind, pbest_fit, pbest_bio, pbest_prod, gbest, gbest_ind, gbest_fit, gbest_bio, gbest_prod] = update_pbest_gbest(particles, indexSelected, fitness, bioValue, ProdValue, pbest, pbest_ind, pbest_fit, pbest_bio, pbest_prod, gbest, gbest_ind, gbest_fit, gbest_bio, gbest_prod)
    % Update pbest and gbest
    for i = 1:size(particles, 1)
        if fitness(i) > pbest_fit(i)
            pbest(i, :) = particles(i, :);
            pbest_ind(i, :) = indexSelected(i, :);
            pbest_fit(i) = fitness(i);
            pbest_bio(i) = bioValue(i);
            pbest_prod(i) = ProdValue(i);
        end
    end
    [max_fit, ind] = max(fitness);
    if max_fit > gbest_fit
        gbest = particles(ind, :);
        gbest_ind = indexSelected(ind, :);
        gbest_fit = max_fit;
        gbest_bio = bioValue(ind);
        gbest_prod = ProdValue(ind);
    end
end

function nbest = CalculateNBest(particles, particleNum, pbest_fitness, neighborhoodSize)
    % CalculateNBest function calculates the neighbor best for each particle.
    %
    % particleNum: The number of particles in the swarm.
    % pbest_fitness: A 1D array containing the fitness values of the pbest of each particle.
    % neighborhoodSize: The number of neighbors considered for each particle.
    %
    % The function returns a 1D array of indices, where the i-th element corresponds to the
    % index of the nbest particle for the i-th particle in the swarm.

    % Preallocate the nbest array
    nbest = zeros(size(particles));

    % Half of the neighborhood size
    halfNeighborhood = floor(neighborhoodSize / 2);

    for i = 1:particleNum
        % Determine the neighborhood indices using modulo operation to handle circular neighborhood
        neighbors = mod((i - 1) + (-halfNeighborhood : halfNeighborhood), particleNum) + 1;

        % Find the index of the neighbor with the best fitness
        [~, maxIdx] = max(pbest_fitness(neighbors));

        % Assign the index of the best neighbor to the nbest array
        nbest(i,:) = particles(maxIdx,:);
    end
end

function mbest = Get_mbest(pbest)
    [l, d] = size(pbest);
    mbest = zeros(1, d);
    
    sums = sum(pbest); % 对pbest进行列求和
    avg = sums / l;
    
    for j = 1:d % 个体维数
        if avg(j) > 0.5
            mbest(j) = 1;
        elseif avg(j) == 0.5
            mbest(j) = rand < 0.5;
        end
    end
end


function Pi = Get_P2(pbesti, gbest, nbesti)
	% Uniform crossover on pbest, gbest and nbest to obtain P
    d = size(gbest,2);
    Pi = zeros(1,d);
    r = rand;
    for i = 1:d
        if r < 1/3
            Pi(i) = pbesti(i);
        elseif r < 2/3
            Pi(i) = gbest(i);
        else
            Pi(i) = nbesti(i);
        end
    end
end

function Pi = Get_P(pbesti, gbest)
    % Uniform crossover on pbest and gbest to obtain P
    d = size(gbest,2);
    Pi = zeros(1,d);
    for i = 1:d
        if rand < 0.5
            Pi(i) = pbesti(i);
        else
            Pi(i) = gbest(i);
        end
    end
end

function Xid = Transf(Pid,Pr)
    dvl =size(Pid,2);
    for i = 1:dvl
        if rand < Pr
            if Pid(i)==0
                Pid(i)=1;
            else
                Pid(i)=0;
            end
        end
    end
    Xid = Pid;
end

