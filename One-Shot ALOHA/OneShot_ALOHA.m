%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization of the simulation
% clear all; % Initialise toutes les variables
close all; % Ferme toutes les fenetres ouvertes
clc; % Clear command window

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ANALYTICAL PARAMETRES
G = 0:0.1:10; % Average attempted RA transmission
% For N = 3
N1 = 3; % Number of RAOs
M1 = 1:1:(10*N1); % Number of devices
S1_number = zeros(30, 1); % Number of successful devices
C1_number = zeros(30, 1); % Number of collided devices
S1_error = zeros(30, 1); % Error of number of successful devices
C1_error = zeros(30, 1); % Error of number of collided devices
% For N = 14
N2 = 14; % Number of RAOs
M2 = 1:1:(10*N2); % Number of devices
S2_number = zeros(140, 1); % Number of successful devices
C2_number = zeros(140, 1); % Number of collided devices
S2_error = zeros(140, 1); % Error of number of successful devices 
C2_error = zeros(140, 1); % Error of number of collided devices
I_MAX = 10; % Number of retransmission
P_S = zeros(41, 1); % Successful probability
T_a = zeros(41, 1); % Mean access delay
P_C = zeros(41, 1); % Collided probability

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 1 : ANALYTICAL AND APPROXIMATION RESULTS
for m = 1:length(M1)
    fprintf("N = 3, M = %d\n", m)
    result_s = 0;
    result_c = 0;
    for k = 0:min(N1, floor(m/2))
        % fprintf("------------------------\n")
        % fprintf("k = %f, M = %f, N = %f\n", k, m, N1)
        % result_s = result_s + (p_k_s('S', k, m, N1) / N1);
        result_s = result_s + (p_k("recursive", 'S', k, m, N1) / N1);
    end
    for k = 1:min(N1, floor(m/2))
        % fprintf("------------------------\n")
        % fprintf("k = %f, M = %f, N = %f\n", k, m, N1)
        % result_c = result_c + (k * p_k_c('C', k, m, N1) / N1);
        result_c = result_c + (k * p_k("recursive", 'C', k, m, N1) / N1);
    end
    S1_number(m) = result_s;
    C1_number(m) = result_c;
end
for m = 1:length(M2)
    fprintf("N = 14, M = %d\n", m)
    result_s = 0;
    result_c = 0;
    for k = 0:min(N2, floor(m/2))
        % fprintf("S, k = %d\n", k)
        result_s = result_s + (p_k("iterative", 'S', k, m, N2) / N2);
    end
    for k = 1:min(N2, floor(m/2))
        % fprintf("C, k = %d\n", k)
        result_c = result_c + (k * p_k("iterative", 'C', k, m, N2) / N2);
    end
    S2_number(m) = result_s;
    C2_number(m) = result_c;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 2 : ABSOLUTE APPROXIMATION ERROR
for m = 1:length(M1)
    S1_error(m) = abs(S1_number(m) - ((m/N1)*exp(-m/N1))) / S1_number(m) * 100;
    C1_error(m) = abs(C1_number(m) - (1-1*exp(-m/N1)-(m/N1)*exp(-m/N1))) / C1_number(m) * 100;
end
for m = 1:length(M2)
    S2_error(m) = abs(S2_number(m) - ((m/N2)*exp(-m/N2))) / S2_number(m) * 100;
    C2_error(m) = abs(C2_number(m) - (1-1*exp(-m/N2)-(m/N2)*exp(-m/N2))) / C2_number(m) * 100;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 3 : ACCESS SUCCESS PROBABILITY
for n = 1:41 % N = 5 to 45
    K = 100; % K = M = 100
    results = 0;
    for i = 1:I_MAX
        results = results + K*exp(-K/(n+4));
        K = K - K*exp(-K/(n+4));
    end
    P_S(n) = results / 100;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 4 : MEAN ACCESS DELAY
for n = 1:41 % N = 5 to 45
    K = 100; % K = M = 100
    results_1 = 0;
    results_2 = 0;
    for i = 1:I_MAX
        results_1 = results_1 + i * K*exp(-K/(n+4));
        results_2 = results_2 + K*exp(-K/(n+4));
        K = K - K*exp(-K/(n+4));
    end
    T_a(n) = results_1 / results_2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 5 : COLLISION PROBABILITY
for n = 1:41 % N = 5 to 45
    K = 100; % K = M = 100
    results_1 = 0;
    results_2 = 0;
    for i = 1:I_MAX
        results_1 = results_1 + ((n+4)-(n+4)*exp(-K/(n+4))-K*exp(-K/(n+4)));
        results_2 = results_2 + (n+4);
        K = K - K*exp(-K/(n+4));
    end
    P_C(n) = results_1 / results_2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SIMULATION PARAMETRES
% S = 1e3; % Number of simulation times
% T = 1e4; % Total simulation time
% C = 5:1:45; % % Number of RAOs
% D = 100; % Number of devices
% I_max = 10; % Number of retransmission
% num_retrans = zeros(1, D); % Number of retransmission
% success_rao = zeros(length(C), length(G)); % Number of successful transmission
% collision_rao = zeros(length(C), length(G)); % Number of collided transmission
% empty_rao = zeros(length(C), length(G)); % Number of idle transmission
% P_G = zeros(T, length(G)); % Proba frame of each G

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATE TOTAL SUCCESSFUL / FAILED TRANSMISSION OF EVERY G
% for s = 1:S % To avoid the error of throughput
%     fprintf('simulation %f time\n', s);
%     for i = 1:length(C) % 2 times
%         P_C = zeros(T, C(i), length(G)); % Proba frame of each slot of each G of each C, 1000*C*101
%         fprintf('Total channel = %f\n', C(i));
%         for j = 1:length(G) % Iterate elements of G
%             % fprintf('frm = %f\n', G(i));
%             P_G(:, j) = poissrnd(G(j), T, 1);
%             for t = 1:T % Iterate the slot of each G
%                 for n = 1:P_G(t, j) % Distribute each frm of G to random C
%                     % randi(i) generate random num between 1 and i
%                     P_C(t, randi(C(i)), j) = P_C(t, randi(C(i)), j) + 1;
%                 end
%                 for c = 1:1:C(i) % C times
%                     if P_C(t, c, j) == 1 % Only 1 frame transmission
%                         success_rao(i, j) = success_rao(i, j) + 1;
%                     elseif P_C(t, c, j) > 1 % More than 1 frame transmission
%                         collision_rao(i, j) = collision_rao(i, j) + 1;
%                     elseif P_C(t, c, j) < 1 % Less than 1 frame transmission
%                         empty_rao(i, j) = empty_rao(i, j) + 1;
%                     end
%                 end
%             end
%         end
%     end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATE AVERAGE NUM OF SUCCESSFUL / COLLISION TRANSMISSION PER T
% average_throughputs = (success_rao / T) / S;
% average_collisions = (collision_rao / T) / S;
% average_empties = (empty_rao / T) / S;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT THE FIGURE
% figure(1);
% title('Analytical and approximation results of N_{S,1}/N and N_{C,1}/N');
% hold on; % Can keep former plotting and add new plotting
% grid on; % Display the grid
% plot(M1/N1, S1_number, 'r:o'); % Result of analysis
% plot(M1/N1, C1_number, 'b:o'); % Result of analysis
% plot(M2/N2, S2_number, 'm:o'); % Result of analysis
% plot(M2/N2, C2_number, 'c:o'); % Result of analysis
% plot(G, G.*exp(-G), 'k-.'); % Result of analysis
% plot(G, 1-exp(-G)-(G).*exp(-G), 'g-.'); % Result of analysis
% xlabel('M/N');
% ylabel('RAOs/N');
% xlim([0 10])
% ylim([0 1])
% legend('N=3 N_{S,1}/N Analytical Model', 'N=3 N_{C,1}/N Analytical Model', 'N=14 N_{S,1}/N Analytical Model', 'N=14 N_{C,1}/N Analytical Model', 'N_{S,1}/N Derived Performance Metric', 'N_{C,1}/N Derived Performance Metric');
% saveas(figure(1), 'figure1.jpg');

% figure(2);
% title('Absolute approximation error of N_{S,1}/N and N_{C,1}/N');
% % hold on; % Can keep former plotting and add new plotting
% grid on; % Display the grid
% semilogy(M1/N1, S1_error, 'r:o', M1/N1, C1_error, 'b:o', M2/N2, S2_error, 'm:o', M2/N2, C2_error, 'c:o'); % Result of analysis
% xlabel('M/N');
% ylabel('Approximation Error (%)');
% xlim([0 10])
% legend('N=3 N_{S,1}/N', 'N=3 N_{C,1}/N', 'N=14 N_{S,1}/N', 'N=14 N_{C,1}/N');
% saveas(figure(2), 'figure2.jpg');

figure(3);
title('Access success probability and its approximation error');
hold on; % Can keep former plotting and add new plotting
grid on; % Display the grid
plot(5:1:45, P_S, 'r-'); % Result of analysis
% plot(, , 'b:o'); % Result of simulation
xlabel('N');
ylabel('Access Success Probability');
xlim([5 45])
ylim([0 1])
legend('Derived Performance Metric');
saveas(figure(3), 'figure3.jpg');

% figure(4);
% title('Mean access delay and its approximation error');
% hold on; % Can keep former plotting and add new plotting
% grid on; % Display the grid
% plot(5:1:45, T_a, 'r-'); % Result of analysis
% % plot(, , 'b:o'); % Result of simulation
% xlabel('N');
% ylabel('Mean Access Delay(unit : Access Cycle)');
% xlim([5 45])
% ylim([0 10])
% legend('Derived Performance Metric');
% saveas(figure(4), 'figure4.jpg');

% figure(5);
% title('Collision probability and its approximation error');
% hold on; % Can keep former plotting and add new plotting
% grid on; % Display the grid
% plot(5:1:45, P_C, 'r-'); % Result of analysis
% % plot(, , 'b:o'); % Result of simulation
% xlabel('N');
% ylabel('Collision Probability');
% xlim([5 45])
% ylim([0 1])
% legend('Derived Performance Metric');
% saveas(figure(5), 'figure5.jpg');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%