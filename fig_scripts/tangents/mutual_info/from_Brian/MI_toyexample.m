%%
% Toy example to demonstrate Mutual Information calculation
% Define input distributions (A, B, C) and corresponding outputs (10 bins)

% Probability of each input category (equal probability)
P_input = [1/3, 1/3, 1/3];

% Output probability distributions conditioned on each input
% (each row corresponds to an input category, columns are output bins)
P_output_given_input = [
    % Bins:  1    2    3    4    5    6    7    8    9    10
             0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1; % A: Uniform
             0.4, 0.2, 0.1, 0.1, 0.05,0.05,0.04,0.03,0.02,0.01; % B: Skewed left
             0.01,0.02,0.03,0.04,0.05,0.05,0.1, 0.1, 0.2, 0.4;  % C: Skewed right
];

%-----------------------------------

%Less MI
%P_output_given_input(2,:) = P_output_given_input(1,:);

%More MI
%P_output_given_input(2,:) = [1 0 0 0 0 0 0 0 0 0];

%-----------------------------------

% Compute joint distribution P(input,output)
P_joint = zeros(3,10);
for i=1:3
    P_joint(i,:) = P_input(i) * P_output_given_input(i,:);
end

% Compute marginal output distribution P(output)
P_output = sum(P_joint,1);

% Calculate Mutual Information (I(X;Y))
I_xy = 0; I_xy2 = 0;
for i=1:3
    for j=1:10
        if P_joint(i,j) > 0
            I_xy = I_xy + P_joint(i,j) * log2(P_joint(i,j)/(P_input(i)*P_output(j)));
            I_xy2 = I_xy2 + P_input(i) * P_output_given_input(i,j)* log2(P_output_given_input(i,j)/P_output(j));
        end
    end
end

% Display results
disp('Joint Probability Distribution P(input,output):');
disp(P_joint);

disp('Marginal Output Probability P(output):');
disp(P_output);

fprintf('Mutual Information (I(X;Y)): %.4f bits\n', I_xy);
fprintf('Mutual Information (I2(X;Y)): %.4f bits\n', I_xy2);
