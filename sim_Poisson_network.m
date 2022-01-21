%% Model parameters
a = 0.05; % packet arrival probability per transmitter
p = 0.50; % medium access probability by the transmitters

lam   = 1.00; % density of users per unit of area
alpha = 5.00; % path-loss exponent
theta = 2.00; % threshold for succesful coomunication (SINR > theta)
W     = 0.00; % average noise power
R     = 0.50; % mean link distance

%% Simulation parameters
STEPS = 1000; % number of simulation time slots

Ar = 1200/lam; % area of the observable region (square)
L  = sqrt(Ar); % side length of Ar

n  = poissrnd(lam*Ar); % number of transmitters inside the region
NQ = n; % number of buffers or queues (one for each transmitter)

x  = L*rand(1,n); % x position of each transmitter in Ar
y  = L*rand(1,n); % y position of each transmitter in Ar

Q   = geornd(1-0.0,1,NQ); % randomize initial conditions of the buffers
r   = raylrnd(sqrt(2/pi)*R,1,NQ); % randomize the link distances (Rayleigh)
% r   = 1.*ones(1,NQ); % set a unique link distance for all links

%% Simulation
rho = zeros(1,STEPS); % initialize variable to measure queue load
p_s = zeros(1,STEPS); % initialize variable to measure success prob.
thr = zeros(1,STEPS); % initialize variable to measure throughput.
for k = 1:STEPS
    Q      = Q + binornd(1,a,1,n); % receive packets
    rho(k) = mean( Q > 0 ); % measure state
    
    TQ    = Q & binornd(1,p,1,n); % transmitting queues
    SINR  = exprnd(1,1,n).*(r.^-alpha)./(W+Interf(n,TQ,x,y,L,alpha));
    Q     = Q - (TQ&(SINR>theta)); % transmit packets
    thr(k)= sum(TQ&(SINR>theta))/Ar;
    p_s(k)= mean( SINR>theta );
end

avg_rho = mean(rho(STEPS/2:STEPS)); % discards the transient
avg_ps  = mean(p_s(STEPS/2:STEPS));
avg_thr = mean(thr(STEPS/2:STEPS));

%% Plots
figure(2)
hold on
plot(1:STEPS,rho,1:STEPS,p_s,1:STEPS,thr)
legend('rho','p_s','thr')

%% Info
fprintf('\nR/L = %f',R/L);
fprintf('\nStationary p_s = %f', avg_ps);
fprintf('\nStationary thr = %f', avg_thr);
fprintf('\nStationary rho = %f\n', avg_rho);

%% The following function calculates the interference used in the SINR
function I = Interf(n,TQ,x,y,L,alpha)
dx = abs(repmat(x',1,n)-repmat(x,n,1));
dx = min( dx, L-dx );
dy = abs(repmat(y',1,n)-repmat(y,n,1));
dy = min( dy, L-dy );
ds = sqrt(dx.^2+dy.^2); % euclidean distance
ds(eye(n)~=0) = Inf; % disconsider the diagonal
I  = sum(repmat(TQ,n,1).*exprnd(1,n,n).*(ds.^-alpha),2)';
end
