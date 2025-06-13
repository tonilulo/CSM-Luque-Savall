clc; clear;

% Data that is 100% dimensionally consistent
p_u = [1 2 3; 4 5 6];         % 2x3
n_u = [0 0 0 100; 0 0 0 200]; % 2x4

p_l = [7 8 9; 10 11 12; 13 14 15];   % 3x3
n_l = [0 0 0 300; 0 0 0 400; 0 0 0 500]; % 3x4

% Define x_u and x_l to match number of rows in p_u and p_l
x_u = zeros(size(p_u,1),1);
x_l = zeros(size(p_l,1),1);

% === Loop version ===
k = 1;
Pe_loop = zeros(3*(length(x_u)+length(x_l)), 3);
for i = 1:length(x_u)
    for j = 1:3
        Pe_loop(k,1) = p_u(i,j);
        Pe_loop(k,2) = n_u(i,4);
        Pe_loop(k,3) = j;
        k = k + 1;
    end
end
for i = 1:length(x_l)
    for j = 1:3
        Pe_loop(k,1) = p_l(i,j);
        Pe_loop(k,2) = n_l(i,4);
        Pe_loop(k,3) = j;
        k = k + 1;
    end
end

% === Vectorized version ===
n1 = length(x_u);
n2 = length(x_l);
Pe_vec = [
    reshape(p_u.', [], 1), repelem(n_u(:,4), 3), repmat((1:3).', n1, 1);
    reshape(p_l.', [], 1), repelem(n_l(:,4), 3), repmat((1:3).', n2, 1)
];

% === Comparison ===
disp("Are Pe_loop and Pe_vec equal?");
disp(isequal(Pe_loop, Pe_vec));  % Should print 1
