%% Drive code
clear; clc;
global M R POPULATION_SIZE GENES;
% M = 1, 2, ... , 18
M = 18
% R = 1, 2, 3
R = 3
% population size 100
POPULATION_SIZE = 100;
% create GENES
GENES = [0,1];

%create initial population
population = zeros(100,54);
for i = 1:100
	gnome = create_gnome();
	population(i,:) = gnome;
end

% iterative step
max_iter = 150;
generation = 0;
while(generation < max_iter)
    % sort the fittness
    fitness_set = zeros(POPULATION_SIZE, 1);
    for i = 1:length(population)
        fitness_set(i) = cal_fitness(population(i,:));
    end
    new_set = sort(fitness_set);
    
    % found or not
    if(cal_fitness(population(1,:)) <= 0)
       found = 1;
    end
   
    % new generation
    new_generation = zeros(100,54);
    
    % 10% of fittest generation goes to the nest generation
    s1 = 10*POPULATION_SIZE/100;
    for i = 1:s1
        new_generation(i,:) = population(find(fitness_set == new_set(i), 1), :);
    end
    
    % from 50% of fittest generation, individuals mate to produce offspring
    s2 = 90*POPULATION_SIZE/100;
    for i = 1:s2
        r = randi([1 50*POPULATION_SIZE/100], 1, 1);
        parent1 = population(find(fitness_set == new_set(r), 1), :);
        r = randi([1 50*POPULATION_SIZE/100], 1, 1);
        parent2 = population(find(fitness_set == new_set(r), 1), :);
        offspring = mate(parent1, parent2);
        new_generation(i+s1, :) = offspring;
    end
    
    population = new_generation;
    generation = generation + 1
    new_set(1)
end

%% related functions

% create chromosome or string of genes
function gnome = create_gnome()
    global GENES;
	gnome = zeros(1, 54);
	for i = 1:54 %M*R = 54
		gnome(i) = mutated_genes();
	end
end

% create random genes for mutation
function gene = mutated_genes()
    global GENES;
	gene = randi(GENES, 1, 1); % randomly chose GENES
end

% perform mating and produce new offspring
function child_chromosome = mate(par1, par2)
    global M R;
	child_chromosome = zeros(1, M*R);
    flag = zeros(1, M);
	for i = 1:R
        for j = 1:M
            p = rand;
            if p < 0.47
                if flag(j) ~= 2 && i ~= R
                    child_chromosome((i-1)*M+j) = par1((i-1)*M+j);
                    if child_chromosome((i-1)*M+j) == 1
                        flag(j) = 2;
                    end
                elseif flag(j) ~= 2 && i == R
                        child_chromosome((i-1)*M+j) = 1;
                        flag(j) = 2;
                end
            elseif p < 0.94
                if flag(j) ~= 2 && i ~= R
                    child_chromosome((i-1)*M+j) = par2((i-1)*M+j);
                    if child_chromosome((i-1)*M+j) == 1
                        flag(j) = 2;
                    end
                elseif flag(j) ~= 2 && i == R
                        child_chromosome((i-1)*M+j) = 1;
                        flag(j) = 2;
                end
            else
                if flag(j) ~= 2 && i ~= R
                    child_chromosome((i-1)*M+j) = mutated_genes();
                    if child_chromosome((i-1)*M+j) == 1
                        flag(j) = 2;
                    end
                elseif flag(j) ~= 2 && i == R
                        child_chromosome((i-1)*M+j) = 1;
                        flag(j) = 2;
                end
            end
        end
	end
end

% claculate fitness score
function fitness = cal_fitness(chromosome)
    global R;
	% saprate minC in 3 parts
	fit_part1 = 0; fit_part2 = 0; fit_rart3 = 0;
    len = length(chromosome);
    
    % parameters
    C_L = [5125330 3379900 2427980; 5125330 3379900 2427980; 5125330 3379900 2427980];
    C_P = [128000 69000 85000; 128000 69000 85000; 128000 69000 85000];
    C_S = [5740000 5740000 5740000; 5075000 5075000 5075000; 3920000 3920000 3920000];
    C_0 = [12400000 12400000 12400000; 10500000 10500000 10500000; 7000000 7000000 7000000];
    % the parameter C_F need to be fixed
    C_F = [1 1 1; 1 1 1; 1 1 1];
    T = 1;
    
    % claculate part1 and part2
	for i = 1:len
        if i > 0 && i <= len/R
            % determine the model of ship
            if i > 0 && i <= 4
                j = 1;
            elseif i > 4 && i <= 12
                j = 2;
            else
                j = 3;
            end
            fit_part1 = fit_part1 + 350/(T/24)*(C_F(j,1)+C_P(j,1)+C_S(j,1))*chromosome(i);
            fit_part2 = fit_part2 + C_L(j,1)*chromosome(i);
        elseif i > len/R && i <= 2*len/R
            % determine the model of ship
            if i-len/R > 0 && i-len/R <= 4
                j = 1;
            elseif i-len/R > 4 && i-len/R <= 12
                j = 2;
            else
                j = 3;
            end
            fit_part1 = fit_part1 + 350/(T/24)*(C_F(j,2)+C_P(j,2)+C_S(j,2))*chromosome(i);
            fit_part2 = fit_part2 + C_L(j,2)*chromosome(i);
        else
            % determine the model of ship
            if i-2*len/R > 0 && i-2*len/R <= 4
                j = 1;
            elseif i-2*len/R > 4 && i-2*len/R <= 12
                j = 2;
            else
                j = 3;
            end
            fit_part1 = fit_part1 + 350/(T/24)*(C_F(j,3)+C_P(j,3)+C_S(j,3))*chromosome(i);
            fit_part2 = fit_part2 + C_L(j,3)*chromosome(i);
        end
    end
    
    % calculate part3
    fit_part3 = 0;
    
    % minC = fitness
    fitness = fit_part1 + fit_part2 + fit_part3;
end
