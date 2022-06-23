%% 1D Simulation of spatial savanna model with precipitation
% Explicit Euler in time, 2D trapezoidal rule in space
% the boundaries are absorbing
% NB fire is modelled as the integral of phi of G
tic
%% Set options for plots/movies
close all
fprintf('\n');
DO_MOVIE = 0; % i.e. write movie to avi file for playback, else just plots
SPACE_TIME_PLOT = 1;
FINAL_PLOT = 0;
%% Numerical method parameters
L = 1; % working on [0,L]
N = 200; % N+1 grid points
delta = L/N;  % spatial discretization parameter
h = 0.1; % time discretisation parameter
n = 2000; % number of time steps
tau = (n-1)*h; % simulations time domain is [0,tau]
%% Function definitions
P_fun = @(x, p_0, p_1) p_0 + p_1.*x;
P_fun_piece = @(x, p_0, p_1, slope) max(min((p_0+p_1/2) + slope*(x-0.5), p_0+p_1),p_0);
mu_p = @(mu, mu_s, p) mu + mu_s.*p;
nu_p = @(nu, nu_s, p) nu + nu_s.*p;
alpha_p = @(alpha, alpha_s, p) alpha + alpha_s.*p;
beta_p = @(beta, beta_s, p) beta + beta_s.*p;
J_F_fun = @(x, a, sigma_F) exp( -( (a-x).^2)/(2*sigma_F^2) )/sqrt(2*pi*sigma_F^2);
J_T_fun = @(x, a, sigma_T) exp( -( (a-x).^2 )/(2*sigma_T^2) )/sqrt(2*pi*sigma_T^2);
W_fun = @(x, a, sigma_W) exp( -( (a-x).^2 )/(2*sigma_W^2) )/sqrt(2*pi*sigma_W^2);
omega = @(w_0_fun, w_1, g, t_1_fun, s_1) w_0_fun + (w_1-w_0_fun)./(1 + exp(-(g-t_1_fun)/s_1));
phi = @(f_0_fun, f_1, g, t_2_fun, s_2) f_0_fun + (f_1-f_0_fun)./(1 + exp(-(g-t_2_fun)/s_2));
%% Model parameters (values from PNAS paper)
p_left=0;
p_right=1;
p_0=p_left;
p_1=(p_right-p_left)/L;
mu = 0.1;
mu_s = 0;
nu = 0.05;
nu_s = 0;

alpha_c = 0.2;
alpha_s = 0.8;
beta_c = 1.5;
beta_s = 0.1;

gamma = 0;

w_0_ref = 0.9;
w_1_ref = 0.4;
t_1_ref = 0.4;
s_1 = 0.01; % s_1 standard value is 0.01
f_0_ref = 0.1;
f_1_ref = 0.9;
t_2_ref = 0.4;
s_2 = 0.05;% s_2 standard value is 0.05

disp = 0.01;
sigma_F = disp; % seed dispersal radius forest trees
sigma_T = disp; % seed dispersal radius savanna trees
sigma_W = disp; % fire spread radius
%% Set up the initial distributions of the cover types on the grid
% each row is one time step of the simulation
% solution is a block matrix [LB; SOL; RB] where LB and RB are fixed
% boundary conditions corresponding to a "Dirichlet type boundary"

G0 = 0.01*rand(1,N+1);%phi(0.95, 0.05, 0:delta:L, 0.5, s_2);
S0 = 0.01*rand(1,N+1);%0.1*ones(1,N+1);
T0 = phi(0.95, 0.05, 0:delta:L, 0.1, 0.01);%0.1*ones(1,N+1);
F0 = 0.01*rand(1,N+1);

CR = G0 + S0 + T0 + F0;
G0 = G0./CR;
S0 = S0./CR;
T0 = T0./CR;
F0 = F0./CR;
LB_G = fliplr(G0(1,2:end));
RB_G = fliplr(G0(1,1:end-1));
LB_S = fliplr(S0(1,2:end));
RB_S = fliplr(S0(1,1:end-1));
LB_T = fliplr(T0(1,2:end));
RB_T = fliplr(T0(1,1:end-1));
LB_F = fliplr(F0(1,2:end));
RB_F = fliplr(F0(1,1:end-1));
G = [0*LB_G G0 0*RB_G];
S = [0*LB_S S0 0*RB_S];
T = [ones(size(LB_T)) T0 0*RB_T];
F = [0*LB_F F0 ones(size(RB_F))];
% compute the convolution for E
X = 0:delta:L;
X_L = X-L;
X_L = X_L(1,1:end-1);
X_R = X+L;
X_R = X_R(1,2:end);
E = ones(1,N+1);
temp_normalise = ones(1,N+1);
% Save tempW matrices to avoid computing them again
tempW = ones(N+1,3*N+1);
Trap = ones(1,3*N+1);
Trap(1,3*N+1)=0.5;
Trap(1,1)=0.5;
for i = 1:N+1
    tempW(i,:) =  W_fun([X_L X X_R],(i-1)*delta,sigma_W);
    temp_normalise(1,i) = sum(tempW(i,:))*delta;
    integrand = tempW(i,:).*(G(1,:)+gamma*(T(1,:)+S(1,:)));
    E(1,i) = sum(integrand.*Trap)*delta;
end
C_W = max(temp_normalise);
E(1,:) = E(1,:)/C_W;
%% Compute the birth and mortality matrices as a function of rainfall
% this computes the standard (deterministic) gradient values
%P_grad = P_fun(X,p_0,p_1); % compute the rainfall gradient along the x-axis
P_grad = P_fun_piece(X,p_0,0.8*p_1,0*p_1);
mu_grad = mu_p(mu, mu_s, P_grad);
nu_grad = nu_p(nu, nu_s, P_grad);
alpha_grad = alpha_p(alpha_c, alpha_s, P_grad);
beta_grad = beta_p(beta_c, beta_s, P_grad);
%% preallocaate some temp variables for efficiency
temp1 = ones(1,N+1);
temp2 = ones(1,N+1);
tempF = ones(N+1,3*N+1);
tempT = ones(N+1,3*N+1);
temp_normalise_T = ones(1,N+1);
temp_normalise_F = ones(1,N+1);
%% pre-calculate the 4D convolution matrices
for k = 1:N+1
    tempF(k,:) = J_F_fun([X_L X X_R],(k-1)*delta,sigma_F);
    temp_normalise_F(1,k) = sum(tempF(k,:))*delta;
    tempT(k,:) = J_T_fun([X_L X X_R],(k-1)*delta,sigma_T);
    temp_normalise_T(1,k) = sum(tempT(k,:))*delta;
end
C_F = max(temp_normalise_F);
C_T = max(temp_normalise_T);
%% The numerical scheme
for i = 2:n
    % compute convolutions for this time step
    progressbar(i,n);
    for k = 1:N+1
        integrand1 = tempF(k,:).*F(i-1,:);
        temp1(1,k) = sum(integrand1.*Trap)*delta;
        integrand2 = tempT(k,:).*T(i-1,:);
        temp2(1,k) = sum(integrand2.*Trap)*delta;
    end
    temp1 = temp1/(C_F);
    temp2 = temp2/(C_T);
    G(i,(N+1):2*N+1) = G(i-1,(N+1):2*N+1) + h*( mu_grad.*S(i-1,(N+1):2*N+1) + nu_grad.*T(i-1,(N+1):2*N+1)...
        - alpha_grad.*temp1.*G(i-1,(N+1):2*N+1) + ...
        phi(f_0_ref, f_1_ref, E(i-1,:), t_2_ref, s_2).*F(i-1,(N+1):2*N+1)...
        - beta_grad.*temp2.*G(i-1,(N+1):2*N+1) );
    S(i,(N+1):2*N+1) = S(i-1,(N+1):2*N+1) + h*( -mu_grad.*S(i-1,(N+1):2*N+1) + beta_grad.*temp2.*G(i-1,(N+1):2*N+1)...
        - alpha_grad.*temp1.*S(i-1,(N+1):2*N+1)...
        - omega(w_0_ref, w_1_ref, E(i-1,:), t_1_ref, s_1).*S(i-1,(N+1):2*N+1) );
    T(i,(N+1):2*N+1) = T(i-1,(N+1):2*N+1) + h*( - nu_grad.*T(i-1,(N+1):2*N+1) ...
        + omega(w_0_ref, w_1_ref, E(i-1,:), t_1_ref, s_1).*S(i-1,(N+1):2*N+1)...
        -alpha_grad.*temp1.*T(i-1,(N+1):2*N+1) );
    
    F(i,(N+1):2*N+1) = ones(1,N+1) - S(i,(N+1):2*N+1) - T(i,(N+1):2*N+1) - G(i,(N+1):2*N+1);
    
    G(i,:) = [0*LB_G G(i,(N+1):2*N+1) 0*RB_G];
    S(i,:) = [0*LB_S S(i,(N+1):2*N+1) 0*RB_S];
    T(i,:) = [ones(size(LB_T)) T(i,(N+1):2*N+1) 0*RB_T];
    F(i,:) = [0*LB_F F(i,(N+1):2*N+1) ones(size(RB_F))];
    
    for k = 1:N+1
        integrand = tempW(k,:).*( G(i,:)+gamma*(T(i,:)+S(i,:)) );
        E(i,k) = delta*sum(integrand.*Trap)./C_W;
    end
    % For numerical stability take max with zero and min with one
    G(i,:) = min(max(G(i,:),0),1);
    S(i,:) = min(max(S(i,:),0),1);
    T(i,:) = min(max(T(i,:),0),1);
    F(i,:) = min(max(F(i,:),0),1);
    E(i,:) = min(max(E(i,:),0),1);
end
% The following output is useful when trying to discern whether or not
% a solution is stationary in time
fprintf('\n');
fprintf('The maximum changes on the grid for each variable at the last time step were:\n');
fprintf(['G: ',num2str(max(abs(G(n,:)-G(n-1,:)))),'\n']);
fprintf(['S: ',num2str(max(abs(S(n,:)-S(n-1,:)))),'\n']);
fprintf(['T: ',num2str(max(abs(T(n,:)-T(n-1,:)))),'\n']);
fprintf(['F: ',num2str(max(abs(F(n,:)-F(n-1,:)))),'\n']);
toc
%% Visualise the solution...
% either in a movie...
if DO_MOVIE
    subplot_vis = struct('cdata',[],'colormap',[]);
    v = VideoWriter(['grass_front_formation_L=',num2str(L),'_n = ',num2str(n),'_h= ',num2str(h),'_delta=',num2str(delta),'_fire=',num2str(sigma_W),'_seedsT=',num2str(sigma_T),'_seedsF=',num2str(sigma_F),'.avi']);
    v.FrameRate = 4;
    open(v)
    f = figure;
    for j = 1:100:n
        f.Name = ['Simulation time: t = ', num2str((j-1)*h)];
        ax1 = subplot(2,2,1);
        plot(X,G(j,N+1:2*N+1),'LineWidth',2);
        title('Grass');
        
        ax2 = subplot(2,2,2);
        plot(X,S(j,N+1:2*N+1),'LineWidth',2);
        title('Saplings');
        
        ax3 = subplot(2,2,3);
        plot(X,T(j,N+1:2*N+1),'LineWidth',2);
        title('Savanna Trees');
        
        ax4 = subplot(2,2,4);
        plot(X,F(j,N+1:2*N+1),'LineWidth',2);
        title('Forest Trees');
        
        xlim([ax1 ax2 ax3 ax4],[0 L])
        ylim([ax1 ax2 ax3 ax4],[0 1])
        writeVideo(v,getframe(gcf))
        subplot_vis(j) = getframe(gcf);
    end
    % fig = figure;
    % movie(fig,subplot_vis,3,1)
    close(v)
    % or just plotting the solution
end
%% Space-time plot of dynamics
if SPACE_TIME_PLOT == 1
    fST = figure;
    fST.Name = ['Evolution over time'];
    
    subplot(2,2,1)
    h1 = pcolor(G);
    shading interp
    title('Grass');
    set(h1, 'EdgeColor', 'none');
    ylabel('Time');
    caxis([0,1])
    set(h1, 'EdgeColor', 'none');
    ylabel('Time');
    caxis([0,1])
    xticks([1 floor((N+1)/2) floor((N+1))]);
    xticklabels({num2str(0), num2str(L/2), num2str(L)});
    yticks([1 floor(n/2) floor(n)]);
    yticklabels({num2str(0), num2str(n*h/2), num2str(n*h)});
    xlabel('\Omega');
    set(gca,'linewidth',4);
    set(gca,'FontSize',20);
    %xticklabels({num2str(0), num2str(L/4), num2str(L/2), num2str(3*L/4), num2str(L)});
    %yticks([0 floor(n/4) floor(n/2) floor(3*n/4) floor(n)]);
    %yticklabels({num2str(0), num2str(n*h/4), num2str(n*h/2), num2str(3*n*h/4), num2str(n*h)});
    
    
    subplot(2,2,2)
    h2 = pcolor(S);
    shading interp
    title('Saplings');
    set(h1, 'EdgeColor', 'none');
    caxis([0,1])
    xticks([1 floor((N+1)/2) floor((N+1))]);
    xticklabels({num2str(0), num2str(L/2), num2str(L)});
    yticks([1 floor(n/2) floor(n)]);
    yticklabels({num2str(0), num2str(n*h/2), num2str(n*h)});
    xlabel('\Omega');
    set(gca,'linewidth',4);
    set(gca,'FontSize',20);
    
    subplot(2,2,3)
    h3 = pcolor(T);
    shading interp
    title('Savanna Trees');
    set(h3, 'EdgeColor', 'none');
    caxis([0,1])
    xticks([1 floor((N+1)/2) floor((N+1))]);
    xticklabels({num2str(0), num2str(L/2), num2str(L)});
    yticks([1 floor(n/2) floor(n)]);
    yticklabels({num2str(0), num2str(n*h/2), num2str(n*h)});
    xlabel('\Omega');
    set(gca,'linewidth',4);
    set(gca,'FontSize',20);
    
    subplot(2,2,4)
    h4 = pcolor(F);
    shading interp
    title('Forest Trees');
    set(h4, 'EdgeColor', 'none');
    caxis([0,1])
    set(h1, 'EdgeColor', 'none');
    caxis([0,1])
    xticks([1 floor((N+1)/2) floor((N+1))]);
    xticklabels({num2str(0), num2str(L/2), num2str(L)});
    yticks([1 floor(n/2) floor(n)]);
    yticklabels({num2str(0), num2str(n*h/2), num2str(n*h)});
    xlabel('\Omega');
    set(gca,'linewidth',4);
    set(gca,'FontSize',20);
end
%%
figure(2);
subplot(2,2,1), plot(0:delta:L,G(end,N+1:2*N+1),'r','LineWidth',2);
xlim([0 L])
ylim([0 1])
xlabel('\Omega');
ylabel('Grass Density');
set(gca,'linewidth',1.25);
set(gca,'FontSize',12);
dim = [.24 .615 .3 .3];
str = sprintf('\\bf\\sigma = %.2f', disp);
annotation('textbox',dim,'String',str,'FitBoxToText','on','LineWidth',1.25,'FontSize',11,'Margin',2.5,'HorizontalAlignment','center');

subplot(2,2,2), plot(0:delta:L,S(end,N+1:2*N+1),'r','LineWidth',2);
xlim([0 L])
ylim([0 1])
xlabel('\Omega');
ylabel('Sapling Density');
set(gca,'linewidth',1.25);
set(gca,'FontSize',12);

subplot(2,2,3), plot(0:delta:L,T(end,N+1:2*N+1),'r','LineWidth',2);
xlim([0 L])
ylim([0 1])
xlabel('\Omega');
ylabel('Savanna Density');
set(gca,'linewidth',1.25);
set(gca,'FontSize',12);

subplot(2,2,4), plot(0:delta:L,F(end,N+1:2*N+1),'r','LineWidth',2);
xlim([0 L])
ylim([0 1])
xlabel('\Omega');
ylabel('Forest Density');
set(gca,'linewidth',1.25);
set(gca,'FontSize',12);
hold on;


%% stacked plots
figure(3);
area(0:delta:L,G(end,N+1:2*N+1)+S(end,N+1:2*N+1)+T(end,N+1:2*N+1)+F(end,N+1:2*N+1),'LineWidth',1.5,'FaceColor',[0 0.39 0]);
hold on;
area(0:delta:L,G(end,N+1:2*N+1)+S(end,N+1:2*N+1)+T(end,N+1:2*N+1),'LineWidth',1.5);
area(0:delta:L,G(end,N+1:2*N+1)+S(end,N+1:2*N+1),'LineWidth',1.5);
area(0:delta:L,G(end,N+1:2*N+1),'LineWidth',1.5,'FaceColor',[0.565 0.933 0.565]);
xlim([0 L]);
ylim([0 1]);
xticks([0 1/2 1]);
xticklabels({num2str(0), num2str(L/2), num2str(L)});
yticks([0 1/2 1]);
yticklabels({num2str(0), num2str(1/2), num2str(1)});
xlabel('\Omega');
ylabel('Density');
legend('Forest','Savanna','Saplings','Grass');
set(gca,'linewidth',1.25);
set(gca,'FontSize',18);

%% Plots of averages versus time
figure(4);
time_interval = 0:h:(n-1)*h;
subplot(2,2,1), plot(time_interval, mean(G(:,N+1:2*N+1),2),'r','LineWidth',1);
xlabel('t');
ylabel('< G(\cdot,t) >');
set(gca,'linewidth',1);
set(gca,'FontSize',20);
axis tight;

subplot(2,2,2), plot(time_interval, mean(S(:,N+1:2*N+1),2),'r','LineWidth',1);
xlabel('t');
ylabel('< S(\cdot,t) >');
set(gca,'linewidth',1);
set(gca,'FontSize',20);
axis tight;

subplot(2,2,3), plot(time_interval,mean(T(:,N+1:2*N+1),2),'r','LineWidth',1);
xlabel('t');
ylabel('< T(\cdot,t) >');
set(gca,'linewidth',1);
set(gca,'FontSize',20);
axis tight;

subplot(2,2,4), plot(time_interval, mean(F(:,N+1:2*N+1),2),'r','LineWidth',1);
xlabel('t');
ylabel('< F(\cdot,t) >');
set(gca,'linewidth',1);
set(gca,'FontSize',20);
axis tight;
