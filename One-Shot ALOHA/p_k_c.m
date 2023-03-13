function proba = p_k_c(mode, k, M, N)
    stp1 = nchoosek(N, k) / (N^M);
    % stp2 = i_k_recursive(mode, 1, k, M, N);
    stp2 = i_k_iterative(mode, k, M, N);

    proba = stp1 * stp2;
end

% To make sure the result of my recursive function is same as the original approach method
% function [proba1, proba2] = p_k_c(k, M, N)
%     stp1 = nchoosek(N, k) / (N^M);
%     stp2_1 = i_k_recursive('C', 1, k, M, N);
%     stp2_2 = 0;
% 
%     if k == 1
%         for i_1 = 2:(M-2*(k-1))
%             if M-i_1 > N-k
%                 stp2_2 = stp2_2 + 0;
%             else
%                 stp2_2 = stp2_2 + nchoosek(M, i_1) * nchoosek(N-k, M-i_1) * factorial(M-i_1);
%             end
%         end
%     elseif k == 2
%         for i_1 = 2:(M-2*(k-1))
%             for i_2 = 2:(M-i_1)
%                 if M-i_1-i_2 > N-k
%                     stp2_2 = stp2_2 + 0;
%                 else
%                     stp2_2 = stp2_2 + nchoosek(M, i_1) * nchoosek(M-i_1, i_2) * nchoosek(N-k, M-i_1-i_2) * factorial(M-i_1-i_2);
%                 end
%             end
%         end
%     elseif k == 3
%         for i_1 = 2:(M-2*(k-1))
%             for i_2 = 2:(M-i_1)
%                 for i_3 = 2:(M-i_1-i_2)
%                     if M-i_1-i_2-i_3 > N-k
%                         stp2_2 = stp2_2 + 0;
%                     else
%                         stp2_2 = stp2_2 + nchoosek(M, i_1) * nchoosek(M-i_1, i_2) * nchoosek(M-i_1-i_2, i_3) * nchoosek(N-k, M-i_1-i_2-i_3) * factorial(M-i_1-i_2-i_3);
%                     end
%                 end
%             end
%         end
%     end
% 
%     proba1 = stp1 * stp2_1;
%     proba2 = stp1 * stp2_2;
%     % fprintf("proba1 = %f, proba2 = %f\n", proba1, proba2)
%     if proba1 == proba2
%         fprintf("True\n")
%     else
%         fprintf("False\n")
%     end
% end

% function [proba] = p_k(k, M, N)
%     x = nchoosek(N, k) / (N^M);
%     y = zeros(k, 1);
%     z = 0;
%     for i_1 = 2:(M-2*(k-1))
%         y(1) = y(1) + nchoosek(M, i_1);
%         for i_2 = 2:(M-i_1)
%             y(2) = y(2) + nchoosek(M-i_1, i_2);
%             for i_3 = 2:(M-i_1-i_2)
%                 y(3) = y(3) + nchoosek(M-i_1-i_2, i_3);
%                 if M-i_1-i_2-i_3 > N-k
%                     z = 0;
%                 else
%                     z = nchoosek(N-k, M-i_1-i_2-i_3) / factorial(M-i_1-i_2-i_3);
%                 end
%             end
%         end
%     end
% 
%     proba = x * y(1) * y(2) * y(3) * z;
% end
