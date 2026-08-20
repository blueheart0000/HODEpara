function [best_solution, best_fitness] = SHADE_call(fobj, lb, ub, dim, max_nfe, pop_size, x0)
    % SHADE_call: Callable wrapper for SHADE 1.1 (Tanabe & Fukunaga, 2013)
    %   [Xbest, fbest] = SHADE_call(fobj, lb, ub, dim, max_nfe, pop_size)
    %   [Xbest, fbest] = SHADE_call(fobj, lb, ub, dim, max_nfe, pop_size, x0)
    %
    % Input:
    %   fobj     - function handle to the objective function
    %   lb       - lower bound (scalar or 1xd vector)
    %   ub       - upper bound (scalar or 1xd vector)
    %   dim      - problem dimension
    %   max_nfe  - maximum number of function evaluations
    %   pop_size - population size (default: 100)
    %   x0       - optional initial solution (1xd row vector)
    %
    % Output:
    %   best_solution - best solution found (1xd row vector)
    %   best_fitness  - best fitness value
    %
    % Based on SHADE 1.1 by Ryoji Tanabe and Alex Fukunaga
    % Reference: Tanabe & Fukunaga, IEEE TEC 2013

%     if nargin < 7
%         x0 = [];
%     end
    if nargin < 6
        pop_size = 100;
    end
    if nargin < 5
        max_nfe = 10000 * dim;
    end

    % Expand scalar bounds to vectors
    if isscalar(lb)
        lb = lb * ones(1, dim);
        ub = ub * ones(1, dim);
    end

    % Build bound matrix for boundConstraint
    lu = [lb; ub];

    % SHADE parameters (from Paper 2: Kaveh & Biabani Hamedani, 2024,
    % using Tanabe's SHADE 1.1.1 code: H=5, p=0.2, archive=NP)
    p_best_rate = 0.2;     % p: pbest rate (SHADE 1.1 default)
    arc_rate = 1.0;        % archive = NP (SHADE 1.1 default)
    memory_size = 5;       % H: memory entries (SHADE 1.1 default)

    % Initialize population
%     pop = repmat(lb, pop_size, 1) + rand(pop_size, dim) .* repmat(ub - lb, pop_size, 1);
%     pop = x0;
%     % If an initial solution x0 is provided, replace one individual
%     if ~isempty(x0)
%         pop(1, :) = x0(:).';
%     end
    pop = x0;
    fitness = zeros(pop_size, 1);
    for i = 1:pop_size
        fitness(i) = fobj(pop(i, :));
    end
    nfes = pop_size;

    best_fitness = min(fitness);
    best_solution = pop(fitness == best_fitness, :);
    best_solution = best_solution(1, :); % ensure row vector

    % Memory initialization
    memory_sf = 0.5 * ones(memory_size, 1);
    memory_cr = 0.5 * ones(memory_size, 1);
    memory_pos = 1;

    % Archive initialization (start empty, per original SHADE)
    archive_NP = round(arc_rate * pop_size);
    archive_pop = zeros(0, dim);

    % Main loop
    while nfes < max_nfe
        [~, sorted_index] = sort(fitness);

        % Select memory indices
        mem_rand_index = ceil(memory_size * rand(pop_size, 1));
        mu_sf = memory_sf(mem_rand_index);
        mu_cr = memory_cr(mem_rand_index);

        % Generate crossover rate CR
        cr = normrnd(mu_cr, 0.1);
        cr(mu_cr == -1) = 0;
        cr = min(max(cr, 0), 1);

        % Generate scaling factor F
        sf = mu_sf + 0.1 * tan(pi * (rand(pop_size, 1) - 0.5));
        bad = (sf <= 0);
        while any(bad)
            sf(bad) = mu_sf(bad) + 0.1 * tan(pi * (rand(sum(bad), 1) - 0.5));
            bad = (sf <= 0);
        end
        sf = min(sf, 1);

        % current-to-pbest/1 mutation
        pNP = max(round(p_best_rate * pop_size), 2);
        randindex = ceil(rand(1, pop_size) .* pNP);
        randindex = max(1, randindex);
        pbest = pop(sorted_index(randindex), :);

        % Select r1 from pop, r2 from pop+archive
        r0 = 1:pop_size;
        popAll = [pop; archive_pop];
        [r1, r2] = gnR1R2(pop_size, size(popAll, 1), r0);

        vi = pop + sf(:, ones(1, dim)) .* (pbest - pop + pop(r1, :) - popAll(r2, :));
        vi = boundConstraint(vi, pop, lu);

        % Crossover (binomial)
        mask = rand(pop_size, dim) > cr(:, ones(1, dim));
        rows = (1:pop_size)';
        cols = floor(rand(pop_size, 1) * dim) + 1;
        jrand = sub2ind([pop_size, dim], rows, cols);
        mask(jrand) = false;
        ui = vi;
        ui(mask) = pop(mask);

        % Evaluate offspring
        children_fitness = zeros(pop_size, 1);
        for i = 1:pop_size
            children_fitness(i) = fobj(ui(i, :));
            nfes = nfes + 1;
            if children_fitness(i) < best_fitness
                best_fitness = children_fitness(i);
                best_solution = ui(i, :);
            end
            if nfes >= max_nfe
                break;
            end
        end
        if nfes >= max_nfe
            break;
        end

        % Selection
        diff = fitness - children_fitness;
        I = (diff > 0); % offspring better than parent

        goodCR = cr(I);
        goodF = sf(I);
        diff_val = diff(I);

        % Update archive with replaced parents (random trim, match LSHADE)
        replaced = pop(I, :);
        n_replaced = size(replaced, 1);
        if n_replaced > 0
            archive_pop = [archive_pop; replaced];
            if size(archive_pop, 1) > archive_NP
                keep_idx = randperm(size(archive_pop, 1), archive_NP);
                archive_pop = archive_pop(keep_idx, :);
            end
        end

        % Update population
        fitness(I) = children_fitness(I);
        pop(I, :) = ui(I, :);

        % Update memory with weighted Lehmer mean
        num_success = length(goodCR);
        if num_success > 0
            sum_diff = sum(diff_val);
            if sum_diff > 0
                w = diff_val / sum_diff;
            else
                w = ones(num_success, 1) / num_success;
            end

            % Lehmer mean for F
            memory_sf(memory_pos) = (w' * (goodF .^ 2)) / (w' * goodF);

            % Lehmer mean for CR (terminal -1 if all successful CR are 0)
            if max(goodCR) == 0
                memory_cr(memory_pos) = -1;
            else
                memory_cr(memory_pos) = (w' * (goodCR .^ 2)) / (w' * goodCR);
            end

            memory_pos = memory_pos + 1;
            if memory_pos > memory_size
                memory_pos = 1;
            end
        end
    end
end
