% i_k_iterative('C', 3, 10, 8) for testing
function z = i_k_iterative(mode, k, M, N)
    % Calculate the sum of nchoosek(n, iterative elements of boundary)
    % f = @(n, boundary) sum(factorial(n)./(factorial(n-boundary).*factorial(boundary)));

    z = 0; % Define output variable (IMPORTANT!!)
    itr_sum_series = zeros(k, 1) + 2; % Initialize array with [2 2 ... 2]

    if k == 0 % Successful i_k calculation starts from 0
        if M > N-k % Combination cannot be found
            z = z + 0; % Combination equals to 0
        else
            z = z + nchoosek(N-k, M) * factorial(M) * (M);
        end
    else % Collided i_k calculation starts from 0
        i_n = k; % Calculation of sum of series starts from the last sigma
        while i_n > 0
            % To avoid the duplicated calculation (e.g. calculate [2 2 2] twice)
            if i_n == k % Always calculating from final result of the number iteration of sum of series
                % fprintf("i_n = %d, %d %d %d\n", i_n, itr_sum_series(1), itr_sum_series(2), itr_sum_series(3))
        
                % All combinations based on the iterative number of each sum of series
                c = 1; % The mutiplied-combinations of each iteration
                for i_k = 1:k
                    % fprintf("n = %d, k = %d\n", M-sum_series(itr_sum_series, 1, i_k), itr_sum_series(i_k))
                    % Don't need this if-else because there's no possible that combination cannot be found
                    % if itr_sum_series(i_k) > M-sum_series(itr_sum_series, 1, i_k)  % Combination cannot be found
                    %     c = c * 0; % Combination equals to 0
                    % else
                    %     % Combinations of inner multiple sum of series
                    %     c = c * nchoosek(M-sum_series(itr_sum_series, 1, i_k-1), itr_sum_series(i_k));
                    % end
                    c = c * nchoosek(M-sum_series(itr_sum_series, 1, i_k-1), itr_sum_series(i_k)); % Combinations of the inner sum of series
                end
           
                % Don't need this if-else because there's no possible that combination cannot be found
                % if M-sum_series(itr_sum_series, 1, k) > N-k  % Combination cannot be found
                %     c = c * 0; % Combination equals to 0
                % else
                %     % Functions of the last sum of series
                %     c = c * nchoosek(N-k, M-sum_series(itr_sum_series, 1, k)) * factorial(M-sum_series(itr_sum_series, 1, k));
                % end

                % Combinations of the last sum of series
                % fprintf("n = %d, k = %d\n", N-k, M-sum_series(itr_sum_series, 1, k))
                if M-sum_series(itr_sum_series, 1, k) > N-k % Combination cannot be found
                    z = z + 0; % Combination equals to 0
                else
                    if mode == 'C'
                        c = c * nchoosek(N-k, M-sum_series(itr_sum_series, 1, k)) * factorial(M-sum_series(itr_sum_series, 1, k));
                    elseif mode == 'S'
                        c = c * nchoosek(N-k, M-sum_series(itr_sum_series, 1, k)) * factorial(M-sum_series(itr_sum_series, 1, k)) * (M-sum_series(itr_sum_series, 1, k));
                    end
                    z = z + c; % The summation of those mutiplied-combinations
                end
                % fprintf("z = %d\n", z)
                % fprintf("------------------\n")
            end
    
            % Go to the next (inner) sum of series
            % sigma_lower = itr_sum_series(i_n); % Increase lower boundary
            % if i_n == 1
            %     sigma_upper = M-2*(k-1);
            % else
            %     sigma_upper = M-sum_series(itr_sum_series, 1, i_n-1);
            % end
            % fprintf("lower = %d, upper = %d\n", sigma_lower, sigma_upper)
    
            % When the sum of serie finishing (lower boundary >= upper boundary)
            % fprintf("lower = %d, upper = %d\n", itr_sum_series(i_n), M-2*(k-i_n)-sum_series(itr_sum_series, 1, i_n-1))
            if itr_sum_series(i_n) >= M-2*(k-i_n)-sum_series(itr_sum_series, 1, i_n-1)
                while itr_sum_series(i_n) >= M-2*(k-i_n)-sum_series(itr_sum_series, 1, i_n-1)
                    itr_sum_series(i_n) = 2; % If outer sigma finished, go inner sigma
                    i_n = i_n - 1;
                    if i_n <= 0 % There's no element 0 in an array
                        break
                    end
                end
            % Represent for-loop of the sigma (sum of serie)
            else
                itr_sum_series(i_n) = itr_sum_series(i_n) + 1;
                i_n = k; % Go to the last sum of series
            end
    
            % fprintf("------------------\n")
        end
    end
end

% Cannot use function "symsum", so creating this one function
function result = sum_series(func, lower, upper)
    result = 0;
    for i = lower:upper
        result = result + func(i);
    end
end
