% i_k_recursive('C', 1, 3, 10, 8) for testing
function z = i_k_recursive(mode, k_itr, k, M, N)
    z = 0; % Define output variable (IMPORTANT!!)
    if k_itr == 0 % If the multiple sum of series of i_k starts from 0
        if k_itr == k % If k_itr is the last sum of the series
            if M > N-k % Combination cannot be found
                z = z + 0; % Combination equals to 0
            else
                % Functions of the last sum of series
                z = z + nchoosek(N-k, M) * factorial(M) * (M);
            end
        else
            % Go to the next (inner) sum of series
            z = 0 + i_k_recursive(mode, k_itr+1, k, M, N);
        end
    elseif k_itr == 1 % If the multiple sum of series of i_k starts from 1
        % fprintf("---------\n")
        % fprintf("k_itr = %d\n", k_itr)
        if k_itr == k  % If k_itr is the last sum of the series
            for i_n = 2:(M-2*(k-1)) % boundary of the first sum of series
                % fprintf("i_%d = %d\n", k_itr, i_n)
                if M-i_n > N-k % Combination cannot be found
                    z = z + 0; % Combination equals to 0
                else
                    if mode == 'C'
                        % Functions of the last sum of series
                        z = z + nchoosek(M, i_n) * nchoosek(N-k, M-i_n) * factorial(M-i_n);
                    elseif mode == 'S'
                        % Functions of the last sum of series
                        z = z + nchoosek(M, i_n) * nchoosek(N-k, M-i_n) * factorial(M-i_n) * (M-i_n);
                    end
                    
                end
                % fprintf("z = %d\n", z)
            end
        else
            for i_n = 2:(M-2*(k-1)) % boundary of the first sum of series
                % fprintf("i_%d = %d\n", k_itr, i_n)
                % Go to the next (inner) sum of series
                z = z + nchoosek(M, i_n) * i_k_recursive(mode, k_itr+1, k, M-i_n, N);
                % fprintf("z = %d\n", z)
            end
        end
        % fprintf("1 : z = %d\n", z)
    elseif k_itr <= k % If the multiple sum of series of i_k starts between 2 and k
        % fprintf("---------\n")
        % fprintf("k_itr = %d\n", k_itr)
        for i_n = 2:M % boundary of the rerst of sum of series
            % fprintf("i_%d = %d\n", k_itr, i_n)
            if k_itr == k % If k_itr is the last sum of the series
                if M-i_n > N-k % Combination cannot be found
                    z = z + 0; % Combination equals to 0
                else
                    if mode == 'C'
                        % Functions of the last sum of series
                        z = z + nchoosek(M, i_n) * nchoosek(N-k, M-i_n) * factorial(M-i_n);
                    elseif mode == 'S'
                        % Functions of the last sum of series
                        z = z + nchoosek(M, i_n) * nchoosek(N-k, M-i_n) * factorial(M-i_n) * (M-i_n);
                    end
                end
            else
                if M-i_n < 2 % Combination cannot be found
                    z = z + 0; % Combination equals to 0
                else
                    % Go to the next (inner) sum of series
                    z = z + nchoosek(M, i_n) * i_k_recursive(mode, k_itr+1, k, M-i_n, N);
                end
            end
            % fprintf("z = %d\n", z)
        end
        % fprintf("2 : z = %d\n", z)
        % fprintf("---------\n")
    end
end

% Only record the iterative number (to memorize the number iterated from each sum of series)
% sum_series(zeros(3, 1), 1, 3, 10)
% function [i_k, k_array] = sum_series(k_array, k_itr, k, M)
%     i_k = 0;
%     fprintf("i_k = %d\n", k_itr)
%     if k_itr == 1
%         for i_n = 2:(M-2*(k-1))
%             % fprintf("i_%d = %d\n", k_itr, i_n)
%             k_array(k_itr) = i_n;
%             k_array(k_itr+1) = sum_series(k_array, k_itr+1, k, M-i_n);
%             fprintf("%d %d %d\n", k_array)
%         end
%         % fprintf("1 : z = %d\n", z)
%     elseif k_itr <= k
%         % fprintf("---------\n")
%         fprintf("k_itr = %d\n", k_itr)
%         for i_n = 2:M
%             % fprintf("i_%d = %d\n", k_itr, i_n)
%             if k_itr == k
%                 k_array(k_itr) = i_n;
%             else
%                 k_array(k_itr) = i_n;
%                 k_array(k_itr+1) = sum_series(k_array, k_itr+1, k, M-i_n);
%             end
%             fprintf("%d %d %d\n", k_array)
%         end
%         % fprintf("2 : z = %d\n", z)
%         % fprintf("---------\n")
%     end
%     % disp(i_k)
% end
