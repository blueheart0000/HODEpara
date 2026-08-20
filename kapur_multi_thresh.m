function fitness = kapur_multi_thresh(thresholds_norm, img_hist, L)
    % KAPUR_MULTI_THRESH Kapur maximum entropy multi-threshold segmentation (maximization)
    %   thresholds_norm : [0,1] normalized threshold vector, length = K
    %   img_hist        : grayscale histogram (1 x L)
    %   L               : number of gray levels (typically 256)
    %   fitness         : negative total entropy (L-SHADE minimizes by default)

    K = length(thresholds_norm);
    t = zeros(1, K);
    for i = 1:K
        t(i) = max(1, min(L-2, round(thresholds_norm(i) * (L-2)) + 1));
    end
    t = sort(t);

    % ensure no duplicate thresholds
    t = unique(t);
    while length(t) < K
        t(end+1) = randi([1, L-2]);
        t = sort(unique(t));
    end

    total_pixels = sum(img_hist);
    prob = img_hist / total_pixels;

    boundaries = [0, t, L-1];
    total_entropy = 0;

    for k = 1:(K+1)
        idx_start = boundaries(k) + 1;
        idx_end   = boundaries(k+1) + 1;
        idx_range = idx_start:idx_end;

        omega_k = sum(prob(idx_range));
        if omega_k > 0 && omega_k < 1
            p_norm = prob(idx_range) / omega_k;
            % avoid log(0)
            p_norm(p_norm == 0) = [];
            if ~isempty(p_norm)
                total_entropy = total_entropy - sum(p_norm .* log(p_norm));
            end
        end
    end

    % L-SHADE minimizes, so return negative value
    fitness = -total_entropy;
end
