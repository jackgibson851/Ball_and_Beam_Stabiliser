%% Script description
% Mini project Mechatronic
% Ball beam 3 simulation script
% Authors: Jack Gibson, Jonathan Ott, Andreas Wiedholz
%% Initialization
A = [0 1 0 0; 0.0290 0 7.0062 0; 0 0 0 1; -3.1991 0 0.1088 0];
B = [0; -4.332; 0; 47.7514];
C = [1 0 0 0;  0 0 1 0];
D = [0];
sys = ss(A, B, C, D);
g = 9.81;
r = 0.0127;

R = 1;  % if R< 0 we get a high SSE and ball falls off beam (value >'L'm)

J = 0.02%0.00935;
m = 0.064; 
L= 0.4255;
step_input_value = 0.1;
ica = 0.1; % initial condition angle alpha
icx = L/2;
tf(sys)
Tf1 = tf([-4.332 0 335], [1 0 -0.1378 0 22.42]);
Tf2 = tf([47.75 0 12.47], [1 0 -0.1378 0 22.42]);
%% System stability (Task 2.a)
[p, z] = pzmap(sys)
disp("System has " + num2str(numel(z)) + " zeros");
disp("System has " + num2str(numel(p)) + " poles");
for i = 1:numel(p)
    if real(p(i)) > 0
        disp("Pole " + num2str(i) + " has a positive real part, so the system is unstable");
    end
end
%% Observability and Controllability
Pc = [B A*B A*A*B A*A*A*B];
if rank(Pc) == 4
    disp("System is controllable");
else
    disp("System is not controllable");
end
Po = [C; C*A; C*A*A; C*A*A*A];
if rank(Po) == 4
    disp("System is observable");
else
    disp("System is not observable");
end

%% Full State Feedback 
use_pp = false; % use pole placement

if use_pp == true
    ica = 0.1;
    zeta = 0.8;
    Ts = 1;
    w_n = Ts/zeta * 4;
    %w_d = 2*pi / 7; % damping frequency
    %w_n = w_d / sqrt(zeta^2-1);
    first_poles = -1.5;%-zeta*w_n + w_n * sqrt(1- zeta^2);    %% Linear system = poles @ (-1.5, -1.5, -20, -20)
    second_poles = -20;%-zeta*w_n - w_n * sqrt(1-zeta^2);
    seperation_factor = 1;
    P = [first_poles; first_poles; second_poles; second_poles];
    K = acker(A,B,P);
else
    Q = C'*C; % diagonal matrix relates to the state vector
    Q(1,1) = 62; % displacement of the ball
    Q(2,2) = 250; % velocity of the ball
    Q(3,3) = 1; % angle of the beam
    Q(4,4) = 150; % angular velocity of the beam
    K = lqr(A,B,Q,R);
end
%% Observer Design 

if use_pp == true
    Li = P .* 10;
    zeta = 0.8;
    w_d = 2*pi / 0.5;
    w_n = 10;%w_d / sqrt(1-zeta^2);
    first_poles = -zeta*w_n + w_n * sqrt(1-zeta^2);
    second_poles = -zeta*w_n - w_n * sqrt(1-zeta^2);
    seperation_factor = 1;
    %Li = [first_poles; first_poles*seperation_factor; second_poles*seperation_factor; second_poles];
    Lo = place(A',C',Li)';
else
   obs_q = 0.5;   % obs_q Linear = 0.5
   full_state_poles = eig(A-B*K)
   Lo = place(A', C', full_state_poles .* 10)';
   %Lo =  lqr (A',C',B*B',obs_q)'; % source https://digitalcollection.zhaw.ch/bitstream/11475/3792/2/2010_Buechi_State-Space-Control-LQR-and-Observer.pdf
end