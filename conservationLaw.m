% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                           %
%    conservationLaw                                                        %
%                                                                           %
%                                                                           %
% OUTPUT: Returns the equations of the conservation laws followed by the    %
%    concentrations of the species of a network. The output variables 'N',  %
%    'W_new', 'conservation_laws_LHS', 'T', and 'model' allow the user to   %
%    view the following, respectively:                                      %
%       - Matrix of reaction vectors of the network                         %
%       - Matrix form of the left-hand side of the conservation laws        %
%       - The left-hand side of the conservation laws                       %
%       - The right-hand side of the conservation laws                      %
%       - Complete network with all the species listed in the 'species'     %
%            field of the structure 'model'                                 %
%                                                                           %
% INPUT: model: a structure, representing the CRN (see README.txt for       %
%    details on how to fill out the structure)                              %
%                                                                           %
% Note: Ideas for some parts of the code was motivated by Soranzo and       %
%          Altafini.                                                        %
%                                                                           %
% Reference: Soranzo N, Altafini C (2009) ERNEST: a toolbox for chemical    %
%    reaction network theory. Bioinform 25(21):2853–2854.                   %
%    https://doi.org/10.1093/bioinformatics/btp513                          %
%                                                                           %
% Created: 28 June 2023                                                     %
% Last Modified: 28 June 2023                                               %
%                                                                           %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



function [N, W_new, conservation_laws_LHS, T, model] = conservationLaw(model)

%
% Step 1: Create a list of all species indicated in the reactions
%

% Initialize list of species
model.species = { };

% Get all species from reactants
for i = 1:numel(model.reaction)
    for j = 1:numel(model.reaction(i).reactant)
        model.species{end+1} = model.reaction(i).reactant(j).species;
    end
end

% Get species from products
for i = 1:numel(model.reaction)
    for j = 1:numel(model.reaction(i).product)
        model.species{end+1} = model.reaction(i).product(j).species;
    end
end

% Get only unique species
model.species = unique(model.species);

% Use lowercase letters for concentration of corresponding species
concentration = lower(model.species);

% Convert the concentration letters to symbolic form
species_concentration = sym(concentration, 'real');



%
% Step 2: Form stoichiometric matrix N
%

% Count the number of species
m = numel(model.species);

% Initialize the matrix of reactant complexes
reactant_complexes = [ ];

% Initialize the matrix of product complexes
product_complexes = [ ];

% Initialize the stoichiometric matrix
N = [ ];

% For each reaction in the model
for i = 1:numel(model.reaction)
  
    % Initialize the vector for the reaction's reactant complex
    reactant_complexes(:, end+1) = zeros(m, 1);
    
    % Fill it out with the stoichiometric coefficients of the species in the reactant complex
    for j = 1:numel(model.reaction(i).reactant)
        reactant_complexes(find(strcmp(model.reaction(i).reactant(j).species, model.species), 1), end) = model.reaction(i).reactant(j).stoichiometry;
    end
    
    % Initialize the vector for the reaction's product complex
    product_complexes(:, end+1) = zeros(m, 1);
    
    % Fill it out with the stoichiometric coefficients of the species in the product complex
    for j = 1:numel(model.reaction(i).product)
        product_complexes(find(strcmp(model.reaction(i).product(j).species, model.species), 1), end) = model.reaction(i).product(j).stoichiometry;
    end
    
    % Create a vector for the stoichiometric matrix: Difference between the two previous vectors
    N(:, end+1) = product_complexes(:, end) - reactant_complexes(:, end);
    
    % If the reaction is reversible
    if model.reaction(i).reversible
      
        % Insert a new vector for the reactant complex: make it same as the product complex
        reactant_complexes(:, end+1) = product_complexes(:, end);
        
        % Insert a new vector for the product complex: make it the same as the reactant complex
        product_complexes(:, end+1) = reactant_complexes(:, end-1);
        
        % Insert a new vector in the stoichiometric matrix: make it the additive inverse of the vector formed earlier
        N(:, end+1) = -N(:, end);
    end
end

% Count the total number of reactions
r = size(N, 2);


%
% Step 3: Form matrix W of conservation laws
%

% Get the left kernel of N
%    Notes:
%       1. This is already a valid conservation law. However, we want one that has only 0s and 1s
%       2. Notice that this is the LEFT kernel; the [right] kernel is null(N,'r') (not transposed)
W = null(N', 'r');

% End the code if there are no conservation laws
if isempty(W)
    W_new = [ ];
    conservation_laws_LHS = [ ];
    T = [ ];
    fprintf('%s has no conservation laws. \n\n', model.id)
    return
end



%
% Step 4: Make sure W has only 0s and 1s
%

% Count the number of conservation laws
conservation_laws_count = size(W, 2);

% Check columns with negative signs
cols_neg = any(W < 0);

% Get the nonnegative columns
W_nonneg = find(cols_neg == 0);

% Check that the nonegative columns are just 0s and 1s
only_0_1 = any(W(:, W_nonneg) < 2);

% Start forming the final W using the nonegative columns with only 0s and 1s
W_new = W(:, W_nonneg(find(only_0_1 == 1)));

% Start taking linear combinations of the columns until nonnegative vectors are found
for i = 2:conservation_laws_count
    
    % Generate all possible combinations taken i columns at a time
    combo = nchoosek(1:conservation_laws_count, i);

    % For each combination
    for j = 1:size(combo, 1)

        % Get the corresponding columns of W and add them
        test = sum(W(:, combo(j, :)), 2);

        % Get the 1s entries in the test vector (for Condition 3 below)
        test_1s = find(test == 1);

        % Condition 1: Check if the test vector is nonnegative
        nonneg = any(test < 0) == 0;

        % Condition 2: Check if all elements of the test vector are 0s and 1s only
        all_0s_1s = any(test > 1) == 0;
        
        % Condition 3: Check if an existing conservation law is in the test vector
        for k = 1:size(W_new, 2)

            % Get all 1s in each existing conservation law
            law_1s = find(W_new(:, k) == 1);

            % Don't consider the test vector if the existing is in the test
            if isequal(intersect(law_1s, test_1s), law_1s) == 1
                not_exist = 0;
                break
            else
                not_exist = 1;
            end
        end

        % If all 3 conditions are satisfied
        if nonneg && all_0s_1s && not_exist

            % Add the new conservation law
            W_new(:, end+1) = test;

            % Stop looking if the number of conservation laws is complete
            if size(W_new, 2) == conservation_laws_count
                break
            end
        end
    end

    % Stop looking if the number of conservation laws is complete
    if size(W_new, 2) == conservation_laws_count
        break
    end
end



%
% Step 5: Display the result
%

% Get the variable letter used as species name
species_variable = regexprep(concentration,'[^a-zA-Z\s]','');

% Remove these from the species concentration variables
species_number_as_string = erase(concentration, species_variable);

% For X1, X2, etc. format, convert strings to numbers
species_number_as_number = cellfun(@str2num, species_number_as_string);

% Rearrange the species
[~, species_order] = sort(species_number_as_number);

% Rearrange the species and rows of W_new according to species number
species_concentration = species_concentration(species_order);
W_new = W_new(species_order, :);

% Initialize the vector of index of nonzero entries
first_nnz_entry = zeros(1, conservation_laws_count);

% Get the first nonzero entry of each column
for i = 1:conservation_laws_count
    first_nnz_entry(i) = find(W_new(:, i), 1, 'first');
end

% Sort the indices
[~, column_order] = sort(first_nnz_entry);

% Rearrange the columns such that the first nonzero entries are in order
W_new = W_new(:, column_order);

% Form the left-hand side of the conservation laws
conservation_laws_LHS = W_new'*species_concentration';

% Create vector of symbolic total amounts
T = sym(strcat('T', string(1:conservation_laws_count)), 'real');

% Create cell array of conservation laws, equal signs, and total amounts
conservation_laws = horzcat(string(conservation_laws_LHS), repmat('=', conservation_laws_count, 1), string(T'));

% Display the result
fprintf('The conservation laws for %s are: \n\n', model.id)
for i = 1:conservation_laws_count
    disp(strjoin(conservation_laws(i, :)))
end

end