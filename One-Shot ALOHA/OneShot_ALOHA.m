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
S1_number = zeros(10*N1, 1); % Number of successful devices
C1_number = zeros(10*N1, 1); % Number of collided devices
S1_error = zeros(10*N1, 1); % Error of number of successful devices
C1_error = zeros(10*N1, 1); % Error of number of collided devices
% For N = 14
N2 = 5; % Number of RAOs
M2 = 1:1:(10*N2); % Number of devices
S2_number = zeros(10*N2, 1); % Number of successful devices
C2_number = zeros(10*N2, 1); % Number of collided devices
S2_error = zeros(10*N2, 1); % Error of number of successful devices 
C2_error = zeros(10*N2, 1); % Error of number of collided devices
I_MAX = 10; % Number of retransmission
P_S = zeros(50, 1); % Successful probability
T_a = zeros(50, 1); % Mean access delay
P_C = zeros(50, 1); % Collided probability

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 1 : ANALYTICAL AND APPROXIMATION RESULTS
for m = 1:length(M1)
    % fprintf("N = 3, M = %d\n", m)
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
    fprintf("M1 = %f, S1 = %f, C1 = %f\n", m, S1_number(m), C1_number(m))
end
for m = 1:length(M2)
    % fprintf("N = 14, M = %d\n", m)
    result_s = 0;
    result_c = 0;
    for k = 0:min(N2, floor(m/2))
        fprintf("S, k = %d\n", k)
        result_s = result_s + (p_k("iterative", 'S', k, m, N2) / N2);
    end
    for k = 1:min(N2, floor(m/2))
        fprintf("C, k = %d\n", k)
        result_c = result_c + (k * p_k("iterative", 'C', k, m, N2) / N2);
    end
    S2_number(m) = result_s;
    C2_number(m) = result_c;
    fprintf("M2 = %f, S2 = %f, C2 = %f\n", m, S2_number(m), C2_number(m))
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
% for n = 1:50 % N = 5 to 45
%     K = 100; % K = M = 100
%     results = 0;
%     for i = 1:I_MAX
%         results = results + K*exp(-K/n);
%         K = K - K*exp(-K/n);
%     end
%     P_S(n) = results / 100;
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 4 : MEAN ACCESS DELAY
% for n = 1:50 % N = 5 to 45
%     K = 100; % K = M = 100
%     results_1 = 0;
%     results_2 = 0;
%     for i = 1:I_MAX
%         results_1 = results_1 + i * K*exp(-K/n);
%         results_2 = results_2 + K*exp(-K/n);
%         K = K - K*exp(-K/n);
%     end
%     T_a(n) = results_1 / results_2;
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 5 : COLLISION PROBABILITY
% for n = 1:50 % N = 5 to 45
%     K = 100; % K = M = 100
%     results_1 = 0;
%     results_2 = 0;
%     for i = 1:I_MAX
%         results_1 = results_1 + (n-n*exp(-K/n)-K*exp(-K/n));
%         results_2 = results_2 + n;
%         K = K - K*exp(-K/n);
%     end
%     P_C(n) = results_1 / results_2;
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SIMULATION PARAMETRES
% S = 1e0; % Number of simulation times
% T = 1e7; % Total simulation time
% C = 50; % % Number of RAOs
% D = 100; % Number of devices
% I_max = 10; % Number of retransmission
% success_rao = zeros(C, 1); % Number of successful RAOs selection
% collision_rao = zeros(C, 1); % Number of collided RAOs selection
% average_throughputs = zeros(C, 1);
% average_collisions = zeros(C, 1);
% average_delay = zeros(C, 1); % Access delay of random access process
% avg_throughputs_err = zeros(C, 1);
% avg_collisions_err = zeros(C, 1);
% avg_delay_err = zeros(C, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATE TOTAL SUCCESSFUL / COLLIDED RAOs SELECTION
% for s = 1:S % To avoid the error of throughput
%     fprintf('simulation %f time\n', s)
%     for c = 1:C
%         fprintf('number of RAOs = %f\n', c)
%         delay_rao = 0; % Total delay of each device from every number of RAOs
%         for t = 1:T % Iterate the slot
%             % fprintf('slot = %f\n', t)
%             device_status = zeros(D, 1); % Wether devices selected RAOs successfully or not
%             device_imax = zeros(D, 1); % Time that devices selected RAOs successfully
%             total_collision = 0;
%             for i = 1:I_max % Random access within I_max
%                 % fprintf('Imax = %f\n', i)
%                 device_rao = zeros(D, 1); % Which RAOs these devices selected
%                 for d = 1:D % All the non-RAO devices starts selecting RAOs
%                     if device_status(d) == 0 % If the devices hasn't selected a RAO
%                         device_rao(d) = randi(c); % Devices select RAOs in normal distribution
%                     end
%                 end
% 
%                 for n = 1:c % Check whether devices selected the same RAO or not
%                     % fprintf('RAO = %f\n', n)
%                     if length(device_rao(device_rao == n)) == 1 % To check each RAO selected by one / multiple device(s)
%                         device_status(device_rao == n) = 1; % Device select RAO successfully
%                         device_imax(device_rao == n) = i; % At which Imax device selected RAO successfully
%                     end
%                     rao_selection = 0;
%                     for d = 1:D
%                         if device_status(d) == 0
%                             if device_rao(d) == n
%                                 rao_selection = rao_selection + 1;
%                             end
%                             if rao_selection > 1
%                                 total_collision = total_collision + 1;
%                                 break
%                             end
%                         end
%                     end % fprintf("----------------------\n")
%                 end % fprintf("----------------------\n")
%             end % fprintf("----------------------\n")
%             success_rao(c) = success_rao(c) + length(find(device_status == 1));
%             collision_rao(c) = collision_rao(c) + total_collision;
%             delay_rao = delay_rao + sum(device_imax);
%         end % fprintf("----------------------\n")
%         average_throughputs(c) = success_rao(c) / D;
%         average_collisions(c) = collision_rao(c) / (c * I_max);
%         average_delay(c) = delay_rao / success_rao(c);
% 
%         avg_throughputs_err(c) = abs(P_S(c) - (average_throughputs(c)/T)) / P_S(c) * 100;
%         avg_collisions_err(c) = abs(P_C(c) - (average_collisions(c)/T)) / P_C(c) * 100;
%         avg_delay_err(c) = abs(T_a(c) - average_delay(c)) / T_a(c) * 100;
%     end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT THE FIGURE
figure(1);
title('Analytical and approximation results of N_{S,1}/N and N_{C,1}/N');
hold on; % Can keep former plotting and add new plotting
grid on; % Display the grid
plot(M1/N1, S1_number, 'r-'); % Result of analysis
plot(M1/N1, C1_number, 'b:o'); % Result of analysis
plot(M2/N2, S2_number, 'm-'); % Result of analysis
plot(M2/N2, C2_number, 'c:o'); % Result of analysis
plot(G, G.*exp(-G), 'k-.'); % Result of analysis
plot(G, 1-exp(-G)-(G).*exp(-G), 'g-.'); % Result of analysis
xlabel('M/N');
ylabel('RAOs/N');
xlim([0 10])
ylim([0 1])
legend('N=3 N_{S,1}/N Analytical Model', 'N=3 N_{C,1}/N Analytical Model', 'N=5 N_{S,1}/N Analytical Model', 'N=5 N_{C,1}/N Analytical Model', 'N_{S,1}/N Derived Performance Metric', 'N_{C,1}/N Derived Performance Metric');
saveas(figure(1), 'figure1.jpg');

figure(2);
title('Absolute approximation error of N_{S,1}/N and N_{C,1}/N');
% hold on; % Can keep former plotting and add new plotting
grid on; % Display the grid
semilogy(M1/N1, S1_error, 'r-', M1/N1, C1_error, 'b:o', M2/N2, S2_error, 'm-', M2/N2, C2_error, 'c:o'); % Result of analysis
xlabel('M/N');
ylabel('Approximation Error (%)');
xlim([0 10])
legend('N=3 N_{S,1}/N', 'N=3 N_{C,1}/N', 'N=5 N_{S,1}/N', 'N=5 N_{C,1}/N');
saveas(figure(2), 'figure2.jpg');

% figure(3);
% yyaxis left
% title('Access success probability and its approximation error');
% hold on; % Can keep former plotting and add new plotting
% grid on; % Display the grid
% plot(5:1:45, P_S(5:45), 'r-'); % Result of analysis
% plot(5:1:45, average_throughputs(5:45)/T, 'b:o'); % Result of simulation
% xlabel('N');
% ylabel('Access Success Probability');
% xlim([5 45])
% ylim([0 1])
% yyaxis right
% plot(5:1:45, avg_throughputs_err(5:45), 'g-.');
% ylabel('Approximation Error (%)');
% ylim([0 100])
% legend('Derived Performance Metric', 'Simulation');
% saveas(figure(3), 'figure3.jpg');

% figure(4);
% yyaxis left
% title('Mean access delay and its approximation error');
% hold on; % Can keep former plotting and add new plotting
% grid on; % Display the grid
% plot(5:1:45, T_a(5:45), 'r-'); % Result of analysis
% plot(5:1:45, average_delay(5:45), 'b:o'); % Result of simulation
% xlabel('N');
% ylabel('Mean Access Delay(unit : Access Cycle)');
% xlim([5 45])
% ylim([0 10])
% yyaxis right
% plot(5:1:45, avg_delay_err(5:45), 'g-.');
% ylabel('Approximation Error (%)');
% ylim([0 100])
% legend('Derived Performance Metric', 'Simulation');
% saveas(figure(4), 'figure4.jpg');

% figure(5);
% yyaxis left
% title('Collision probability and its approximation error');
% hold on; % Can keep former plotting and add new plotting
% grid on; % Display the grid
% plot(5:1:45, P_C(5:45), 'r-'); % Result of analysis
% plot(5:1:45, average_collisions(5:45)/T, 'b:o'); % Result of simulation
% xlabel('N');
% ylabel('Collision Probability');
% xlim([5 45])
% ylim([0 1])
% yyaxis right
% plot(5:1:45, avg_collisions_err(5:45), 'g-.');
% ylabel('Approximation Error (%)');
% ylim([0 2])
% legend('Derived Performance Metric', 'Simulation');
% saveas(figure(5), 'figure5.jpg');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%