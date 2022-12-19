%% 1D Simulation of spatial savanna model with precipitation gradient
% Explicit Euler in time, trapezoidal rule in space
tic
%clear all
%close all

%% Set options for plots/movies
%close all
fprintf('\n');
DO_MOVIE = 0; % i.e. write movie to avi file for playback, else just plots
SPACE_TIME_PLOT = 0;
SPATIAL_AVERAGES = 0;
STACKED_PLOTS = 1;
FINAL_PLOT = 0;
PHASE_SPACE_PLOTS = 0;
%% Numerical method parameters
L = 1; % working on [0,L]
N = 500; % N+1 grid points
delta = L/N;  % spatial discretization parameter
h = 0.1; % time discretisation parameter
n = 2000; % number of time steps
tau = (n-1)*h; % simulations time domain is [0,tau]
BC = 'reflecting'; % current options: 'reflecting', 'open', 'periodic'
rng(18502984,'twister');
%% Function definitions
P_fun = @(x, p_0, p_1) p_0 + p_1.*x/L;
P_fun_piece = @(x, p_0, p_1, slope,L) max(min( (p_0+p_1/2) + slope*(x-L/2)/L, p_0+p_1),p_0);
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
p_1=(p_right-p_left);

grad_slope = 1; % modifies the path through the gradient

mu = 0.1;
mu_s = 0;
nu = 0.05;
nu_s = 0;

% alpha_left=1.2;
% alpha_right=1.2;
% 
% beta_left=0.2;
% beta_right=0.2;

alpha_c = 0.8;
alpha_s = 0.5;
beta_c = 0.15;
beta_s = 0.1; %beta_right-beta_left;

gamma = 0;

w_0_ref = 0.9;
w_1_ref = 0.4;
t_1_ref = 0.4;
s_1 = 0.01; % s_1 standard value is 0.01
f_0_ref = 0.1;
f_1_ref = 0.9;
t_2_ref = 0.4;
s_2 = 0.05;% s_2 standard value is 0.05

disp = 0.025;
sigma_F = disp; % seed dispersal radius forest trees
sigma_T = disp; % seed dispersal radius savanna trees
sigma_W = disp; % fire spread radius
%% Set up the initial distributions of the cover types on the grid
% each row is one time step of the simulation
% solution is a block matrix [LB; SOL; RB]

% forest IC
G0 = 0.1*ones(1,N+1);%phi(0.95, 0.05, 0:delta:L, 0.5, s_2);
S0 = 0.05*ones(1,N+1);%0.1*ones(1,N+1);
T0 = 0.05*ones(1,N+1);%phi(0.95, 0.05, 0:delta:L, 0.1, 0.01);%0.1*ones(1,N+1);
F0 = 0.8*ones(1,N+1);

% Savanna with no pinning IC
% G0 = ones(1,N+1)-0.05*(0:delta:L);%phi(0.95, 0.05, 0:delta:L, 0.5, s_2); % 0.5*ones(1,N+1);
% S0 = 0.05*ones(1,N+1);%0.1*ones(1,N+1);
% T0 = phi(0.4, 0.05, 0:delta:L, 0.3, 0.01);%%0.1*ones(1,N+1);
% F0 = zeros(1,N+1)+0.15*(0:delta:L);

% Savanna with pinning IC
% G0 = 0.5*ones(1,N+1);
% S0 = 0.05*ones(1,N+1);
% T0 = 0.1*ones(1,N+1);
% F0 = (0:delta:L);


CR = G0 + S0 + T0 + F0;
G0 = G0./CR;
S0 = S0./CR;
T0 = T0./CR;
F0 = F0./CR;
switch BC
    case 'reflecting'
        LB_G = fliplr(G0(1,2:end));
        RB_G = fliplr(G0(1,1:end-1));
        LB_S = fliplr(S0(1,2:end));
        RB_S = fliplr(S0(1,1:end-1));
        LB_T = fliplr(T0(1,2:end));
        RB_T = fliplr(T0(1,1:end-1));
        LB_F = fliplr(F0(1,2:end));
        RB_F = fliplr(F0(1,1:end-1));
    case 'open'
        LB_G = zeros(1,N);
        RB_G = zeros(1,N);
        LB_S = zeros(1,N);
        RB_S = zeros(1,N);
        LB_T = zeros(1,N);
        RB_T = zeros(1,N);
        LB_F = zeros(1,N);
        RB_F = zeros(1,N);
    case 'periodic'
        LB_G = G0(1,2:end);
        RB_G = G0(1,1:end-1);
        LB_S = S0(1,2:end);
        RB_S = S0(1,1:end-1);
        LB_T = T0(1,2:end);
        RB_T = T0(1,1:end-1);
        LB_F = F0(1,2:end);
        RB_F = F0(1,1:end-1);
    otherwise
        error('boundary condition mispecified...')
end
G = zeros(n,3*N+1);
S = zeros(n,3*N+1);
T = zeros(n,3*N+1);
F = zeros(n,3*N+1);
G(1,:) = [LB_G G0 RB_G];
S(1,:) = [LB_S S0 RB_S];
T(1,:) = [LB_T T0 RB_T];
F(1,:) = [LB_F F0 RB_F];
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
%P_grad = P_fun(X,p_0,p_1); % compute the rainfall gradient along the x-axis
P_grad = phi(0, 1, X, 0.5, 0.1);%P_fun_piece(X,p_0,p_1,grad_slope,L)+0.0*(rand(1,N+1)-0.5);
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
    % now need to update the extended parts of the solution
    switch BC
        case 'reflecting'
            G(i,1:N) = fliplr(G(i,N+2:2*N+1));
            G(i,2*N+2:3*N+1) = fliplr(G(i,(N+1):2*N));
            S(i,1:N) = fliplr(S(i,N+2:2*N+1));
            S(i,2*N+2:3*N+1) = fliplr(S(i,(N+1):2*N));
            T(i,1:N) = fliplr(T(i,N+2:2*N+1));
            T(i,2*N+2:3*N+1) = fliplr(T(i,(N+1):2*N));
            F(i,1:N) = fliplr(F(i,N+2:2*N+1));
            F(i,2*N+2:3*N+1) = fliplr(F(i,(N+1):2*N));
        case 'open'
        case 'periodic'
            G(i,1:N) = G(i,N+2:2*N+1);
            G(i,2*N+2:3*N+1) = G(i,(N+1):2*N);
            S(i,1:N) = S(i,N+2:2*N+1);
            S(i,2*N+2:3*N+1) = S(i,(N+1):2*N);
            T(i,1:N) = T(i,N+2:2*N+1);
            T(i,2*N+2:3*N+1) = T(i,(N+1):2*N);
            F(i,1:N) = F(i,N+2:2*N+1);
            F(i,2*N+2:3*N+1) = F(i,(N+1):2*N);
        otherwise
            error('boundary condition mispecified...')
    end
    
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
%% Gradient plots
figure;
plot(0:delta:L, P_grad,'LineWidth',3);
hold on;
%plot(0:delta:L, beta_grad,'LineWidth',3);
xlim([0 L]);
ylim([-0.25 1.25]);
xlabel('\Omega');
%legend('\alpha(x)','\beta(x)','Location','NorthWest');
set(gca,'linewidth',1.5);
set(gca,'FontSize',30);

%% Phase-space plots of the dynamics at certain grid points (0.25, 0.5 & 0.75)
if PHASE_SPACE_PLOTS == 1
    figure;
    plot3(G(2500:end,ceil(1.25*N+1)),T(2500:end,ceil(1.25*N+1)),F(2500:end,ceil(1.25*N+1)),'LineWidth',3);
    hold on;
    plot3(G(2500:end,ceil(1.5*N+1)),T(2500:end,ceil(1.5*N+1)),F(2500:end,ceil(1.5*N+1)),'LineWidth',3);
    plot3(G(2500:end,ceil(1.75*N+1)),T(2500:end,ceil(1.75*N+1)),F(2500:end,ceil(1.75*N+1)),'LineWidth',3);
    title(['Slope = ',num2str(grad_slope)]);
    legend('x = 0.25','x=0.5','x=0.75');
    xlabel('Grass')
    ylabel('Savanna')
    zlabel('Forest')
    xlim([0 1]), ylim([0 1]), zlim([0 1]);
    grid on;
    set(gca,'linewidth',2);
    set(gca,'FontSize',14);
    view(135,55);
end
%% Space-time plot of dynamics
if SPACE_TIME_PLOT == 1
    % For 6B gradient, subcritical Hopf at P=0.1802, heteroclinic to
    % saddle at P = 0.43328...
    [tempon,osc_onset] = min(abs(P_grad-0.1802));
    [tempoff,osc_offset] = min(abs(P_grad-0.43328));
    
    fST = figure;
    fST.Name = ['Evolution over time'];

    subplot(2,2,1)
    h1 = pcolor(G(:,N+1:2*N+1));
    shading interp
    %title('Grass');
    set(h1, 'EdgeColor', 'none');
    ylabel('Time');
    caxis([0,1])
    xticks([1 floor((N+1)/2) floor((N+1))]);
    xticklabels({num2str(0), num2str(L/2), num2str(L)});
    yticks([1 floor(n/2) floor(n)]);
    yticklabels({num2str(0), num2str(n*h/2), num2str(n*h)});
    xlabel('\Omega');
    hold on;
%     xline(osc_onset,'--r','LineWidth',1);
%     xline(osc_offset,'--r','LineWidth',1);
    set(gca,'linewidth',4);
    set(gca,'FontSize',20);
    %xticklabels({num2str(0), num2str(L/4), num2str(L/2), num2str(3*L/4), num2str(L)});
    %yticks([0 floor(n/4) floor(n/2) floor(3*n/4) floor(n)]);
    %yticklabels({num2str(0), num2str(n*h/4), num2str(n*h/2), num2str(3*n*h/4), num2str(n*h)});
    
    subplot(2,2,2)
    h2 = pcolor(S(:,N+1:2*N+1));
    shading interp
    %title('Saplings');
    set(h1, 'EdgeColor', 'none');
    caxis([0,1])
    xticks([1 floor((N+1)/2) floor((N+1))]);
    xticklabels({num2str(0), num2str(L/2), num2str(L)});
    yticks([1 floor(n/2) floor(n)]);
    yticklabels({num2str(0), num2str(n*h/2), num2str(n*h)});
    xlabel('\Omega');
    hold on;
%     xline(osc_onset,'--r','LineWidth',1);
%     xline(osc_offset,'--r','LineWidth',1);
    set(gca,'linewidth',4);
    set(gca,'FontSize',20);
    
    subplot(2,2,3)
    h3 = pcolor(T(:,N+1:2*N+1));
    shading interp
    %title('Savanna Trees');
    set(h3, 'EdgeColor', 'none');
    caxis([0,1])
    xticks([1 floor((N+1)/2) floor((N+1))]);
    xticklabels({num2str(0), num2str(L/2), num2str(L)});
    yticks([1 floor(n/2) floor(n)]);
    yticklabels({num2str(0), num2str(n*h/2), num2str(n*h)});
    xlabel('\Omega');
    hold on;
%     xline(osc_onset,'--r','LineWidth',1);
%     xline(osc_offset,'--r','LineWidth',1);
    set(gca,'linewidth',4);
    set(gca,'FontSize',20);
    
    subplot(2,2,4)
    h4 = pcolor(F(:,N+1:2*N+1));
    hold on;
    
    shading interp
    %title('Forest Trees');
    set(h4, 'EdgeColor', 'none');
    caxis([0,1])
    set(h1, 'EdgeColor', 'none');
    caxis([0,1])
    xticks([1 floor((N+1)/2) floor((N+1))]);
    xticklabels({num2str(0), num2str(L/2), num2str(L)});
    yticks([1 floor(n/2) floor(n)]);
    yticklabels({num2str(0), num2str(n*h/2), num2str(n*h)});
    xlabel('\Omega');
    hold on;
%     xline(osc_onset,'--r','LineWidth',1);
%     xline(osc_offset,'--r','LineWidth',1);
    set(gca,'linewidth',4);
    set(gca,'FontSize',20);
    
end
%% Plot solution at the final time point
if FINAL_PLOT == 1
    figure;
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
end
%% stacked plots
if STACKED_PLOTS == 1
    figure;
    area(0:delta:L,G(end,N+1:2*N+1)+S(end,N+1:2*N+1)+T(end,N+1:2*N+1)+F(end,N+1:2*N+1),'LineWidth',0.001,'FaceColor',[0 0.39 0]);
    hold on;
    area(0:delta:L,G(end,N+1:2*N+1)+S(end,N+1:2*N+1)+T(end,N+1:2*N+1),'LineWidth',0.01);
    area(0:delta:L,G(end,N+1:2*N+1)+S(end,N+1:2*N+1),'LineWidth',0.001);
    area(0:delta:L,G(end,N+1:2*N+1),'LineWidth',0.001,'FaceColor',[0.565 0.933 0.565]);
    xlim([0 L]);
    ylim([0 1]);
    %xticks([0 1/2 1]);
    %xticklabels({num2str(0), num2str(L/2), num2str(L)});
    %yticks([0 1/2 1]);
    %yticklabels({num2str(0), num2str(1/2), num2str(1)});
    xlabel('\Omega');
    %ylabel('Density');
    legend('Forest','Savanna','Saplings','Grass');
    set(gca,'linewidth',1.5);
    set(gca,'FontSize',30);
end
%% Plots of averages versus time
if SPATIAL_AVERAGES == 1
    figure;
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
end
%%
% figure
% plot(alpha_grad,'.');
% 
% figure
% plot(beta_grad,'.');