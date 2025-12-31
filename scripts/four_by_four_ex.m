
A = [-1 0 -1 -3;
    1 -1 -3 3;
    0 2 1 -5;
    0 -1 4 -2];
B = eye(4);
C = eye(4);
D = zeros(4);

nu = 0; nw = 4; ny = 0; nr = 4;

% ... build or load P here (ss/tf/zpk) ...
idx.u = 1:nu;                % input columns for u
idx.w = nu+(1:nw);           % input columns for w
idx.y = 1:ny;                % output rows for y
idx.r = ny+(1:nr);           % output rows for r

% LTI model
P = ss(A,B,C,D);          % continuous-time by default
P.InputName  = {'w1','w2','w3','w4'};
P.OutputName = {'r1','r2','r3','r4'};

[M, ~] = build_M_from_ss(P, idx);       % <-- this M is G_rw


if false
    % Get the impulse response data
    [y, t] = impulse(P);
    
    % Now you can plot each response individually if needed
    for i = 1:size(y, 3) % Loop over inputs
        figure;
        plot(t, squeeze(y(:,:,i))); % Squeeze to remove singleton dimensions
        title(['Impulse Response for Input ' num2str(i)]);
        xlabel('Time (s)');
        ylabel('Response');
    end
end
