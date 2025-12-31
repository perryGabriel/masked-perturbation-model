% Define the plant P(s)
P_ = tf(3, [3 1]);

% Define the controller K(s)
K = tf(3, [3 10]);

% Calculate the open-loop transfer function
den = (1 + 8*P_*K);
num = [P_ P_*K; -8*P_*K K];

M = num / den;

nu = 0; nw = 2; ny = 0; nr = 2;

% ... build or load P here (ss/tf/zpk) ...
idx.u = 1:nu;                % input columns for u
idx.w = nu+(1:nw);           % input columns for w
idx.y = 1:ny;                % output rows for y
idx.r = ny+(1:nr);           % output rows for r

% plot(impulse(M(2,1)))
wgrid = logspace(-2, 2, 2001);




wmask_1 = [1 0];    % length q
rmask_1 = [1 0];    % length p
[M_d1, Delta_1, wstar_a1, sigmax_a1, info_a1] = process_masks(M, wmask_1, rmask_1, wgrid);

wmask_2 = [0 1];    % length q
rmask_2 = [0 1];    % length p
[M_d2, Delta_2, wstar_a2, sigmax_a2, info_a2] = process_masks(M, wmask_2, rmask_2, wgrid);

wmask_3 = [0 1];    % length q
rmask_3 = [1 0];    % length p
[M_d3, Delta_3, wstar_a3, sigmax_a3, info_a3] = process_masks(M, wmask_3, rmask_3, wgrid);


wmask_1 = [1 1];    % length q
rmask_1 = [1 0];    % length p
[M_1, ~, wstar_d1, sigmax_d1, info_d1] = process_masks(M, rmask_1, wmask_1, wgrid);

wmask_2 = [1 0];    % length q
rmask_2 = [1 1];    % length p
[M_2, ~, wstar_d2, sigmax_d2, info_d2] = process_masks(M, rmask_2, wmask_2, wgrid);

wmask_3 = [0 1];    % length q
rmask_3 = [1 1];    % length p
[M_3, ~, wstar_d3, sigmax_d3, info_d3] = process_masks(M, rmask_3, wmask_3, wgrid);

wmask_4 = [1 0];    % length q
rmask_4 = [1 0];    % length p
[M_4, ~, wstar_d4, sigmax_d4, info_d4] = process_masks(M, rmask_4, wmask_4, wgrid);

% print_info(M_1, Delta_1, wstar_1, sigmax_1, info_1);
% print_info(M_1, Delta_2, wstar_2, sigmax_2, info_2);
% print_info(M_1, Delta_3, wstar_3, sigmax_3, info_3);

Delta_1 = [tf([2.8 -1.62],[1 0.58]) 0; 0 0];
Delta_2 = [0 0; 0 tf([4.7 -0.8],[1 0.001])];
Delta_3 = [0 tf([-1.2 0.7],[1 0.6]); 0 0];

% Create a cell array for plants and perturbations
perturbations = {Delta_1, Delta_2,Delta_3};
plants = {M_1, M_2, M_3, M_4};

% Loop through each combination of plant and perturbation
for i = 1:length(plants)
    for j = 1:length(perturbations)
        % Define the closed-loop system with feedback
        T = feedback(plants{i} * perturbations{j}, eye(2));

        % Check the poles of the closed-loop system
        poles = pole(T);

        % Display results
        fprintf('Closed-loop system with M_%d and Delta_%d:\n', i, j);
        disp('Poles:');
        disp(poles);

        % Check stability
        if all(real(poles) < 0)
            fprintf('Stability: Stable\n\n');
        else
            fprintf('Stability: Unstable\n\n');
        end
    end
end
