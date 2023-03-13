function proba = p_k(methode, mode, k, M, N)
    stp1 = nchoosek(N, k) / (N^M);
    if methode == "iterative"
        if mode == 'C'
            stp2 = i_k_iterative(mode, k, M, N);
        elseif mode == 'S'
            stp2 = i_k_iterative(mode, k, M, N);
        end
    elseif methode == "recursive"
        if mode == 'C'
            stp2 = i_k_recursive(mode, 1, k, M, N);
        elseif mode == 'S'
            stp2 = i_k_recursive(mode, 0, k, M, N);
        end
    end
    proba = stp1 * stp2;
end