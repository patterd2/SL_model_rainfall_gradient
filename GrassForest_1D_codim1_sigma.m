tic
%% Numerical method parameters
L = 1; % working on [0,L]
N = 300; % N+1 grid points
delta = L/N;  % spatial discretization parameter
%% Function definitions
P_fun = @(x, p_0, p_1) p_0 + p_1.*x;
alpha_p = @(alpha, alpha_s, p) alpha + alpha_s.*p;
J_F_fun = @(x, a, sigma_F) exp( -( (a-x).^2)/(2*sigma_F^2) )/(sqrt(2*pi*sigma_F^2)); %1./(pi*sigma_F*(1+((x-a).^2)/sigma_F^2));%
W_fun = @(x, a, sigma_W) exp( -( (a-x).^2 )/(2*sigma_W^2) )/(sqrt(2*pi*sigma_W^2)); %1./(pi*sigma_W*(1+((x-a).^2)/sigma_W^2));%
phi = @(f_0_fun, f_1, g, t_2_fun, s_2) f_0_fun + (f_1-f_0_fun)./(1 + exp(-(g-t_2_fun)/s_2));
%% Model parameters (values from PNAS paper)
p_left=0;
p_right=1;
p_0=p_left;
p_1=(p_right-p_left)/L;

alpha = 0.5;
alpha_s = 1.25;

f_0_ref = 0.1;
f_1_ref = 0.9;
t_2_ref = 0.4;
s_2 = 0.05; % s_2's standard value is 0.05

disp = 0.1:0.01:0.2;
%disp = fliplr(disp);
branches = 1;
sol_norm = zeros(branches,length(disp)); error = zeros(branches,length(disp));
principal_eig = zeros(branches,length(disp));
error_tolerance = 1e-10; % for the nonlinear solver
fig_count = 0;
%G = load('unstable_branch.mat'); % fill in unstable branch of the saddle
%G = G.G;
for jj = 1:branches % 1 for high grass IC, 2 for low grass IC,
    % 4 for mixed spatial stripes between the two species, 3 for a sigmoidal transition
    for count = 1:length(disp)
        fprintf(['progress = ',num2str((count-1)/(length(disp))*100),'%']);
        fprintf('\n');
        fig_count = fig_count + 1;
        sigma_F = disp(count); % seed dispersal radius forest trees
        sigma_W = 0.01; % fire spread radius
        %% Set up the initial distributions of the cover types on the grid
        if jj == 19 && count == 1
            G = 1-0.15*(0:delta:L);%(1+cos(12*pi*(0:delta:L)));
        elseif jj == 71 && count == 1
            G = 0.05*(1+cos(2*pi*(0:delta:L)));
        elseif jj == 11 && count == 1
            G = 0.15*(1 + cos(4*pi*(0:delta:L)) + sin(6*pi*(0:delta:L)) );
        elseif jj == 17 && count == 1
            G = phi(0.95, 0.05, 0:delta:L, 0.5, 0.05); % front pinned branch first
        elseif jj == 3 && count == 1
            G = phi(0.95, 0.6, 0:delta:L, 0.55, 0.15); % trying to find unstable branch connecting saddles
        elseif jj == 1 && count == 1
            G = 0.44*(0:delta:L).^2 - 0.62*(0:delta:L)+0.29;%phi(0.4, 0.1, 0:delta:L, 0.2, 0.1);%0.2 - 0.1*(0:delta:L); %.*(1+cos(12*pi*(0:delta:L)));
        end
        %% plot the initial condition for this simulation
        figure(fig_count);
        subplot(1,2,1), plot(0:delta:L,G,'-.k','LineWidth',2);
        hold on;
        xlim([0 L]); ylim([0 1]);
        set(gca,'linewidth',1.25); set(gca,'FontSize',18);
        xlabel('\Omega'); ylabel('Density');
        set(gca,'linewidth',1.25); set(gca,'FontSize',12); hold on;
        %% set up the spatial grid and pre calculate the convolution matrices
        X = 0:delta:L;
        X_L = X-L;
        X_L = X_L(1,1:end-1);
        X_R = X+L;
        X_R = X_R(1,2:end);
        temp_normalise = ones(1,N+1);
        % Save tempW matrices to avoid computing them again
        tempW = ones(N+1,3*N+1);
        Trap = ones(1,3*N+1);
        Trap(1,3*N+1)=0.5;
        Trap(1,1)=0.5;
        for i = 1:N+1
            tempW(i,:) =  W_fun([X_L X X_R],(i-1)*delta,sigma_W);
            temp_normalise(1,i) = sum(tempW(i,:))*delta;
        end
        C_W = max(temp_normalise);
        % Compute the birth and mortality matrices as a function of rainfall
        P_grad = P_fun(X,p_0,p_1); % compute the rainfall gradient along the x-axis
        alpha_grad = alpha_p(alpha, alpha_s, P_grad);
        tempF = ones(N+1,3*N+1);
        temp_normalise_F = ones(1,N+1);
        %% pre-calculate the convolution matrices
        for k = 1:N+1
            tempF(k,:) = J_F_fun([X_L X X_R],(k-1)*delta,sigma_F);
            temp_normalise_F(1,k) = sum(tempF(k,:))*delta;
        end
        C_F = max(temp_normalise_F);
        %% Solve the nonlinear problem f(x) = 0 for the equilibrium solution
        % first version is the "reflecting" boundary condition
        %f = @(x) - alpha_grad.*((delta*sum( tempF.*(1 - repmat([fliplr(x(1:N)) x fliplr(x(2:N+1))],N+1,1)).*Trap, 2))'/C_F).*x + ...
        %    phi(f_0_ref, f_1_ref, ((delta*sum( tempW.*(repmat([fliplr(x(1:N)) x fliplr(x(2:N+1))],N+1,1)).*Trap, 2))'/C_W), t_2_ref, s_2).*(1-x);
        % second version is the "open" boundary condition
        f = @(x) - alpha_grad.*((delta*sum( tempF.*(repmat([fliplr(0*x(1:N)) 1-x 0*fliplr(x(2:N+1))],N+1,1)).*Trap, 2))'/C_F).*(x) + ...
            phi(f_0_ref, f_1_ref, ((delta*sum( tempW.*((repmat([0*fliplr(x(1:N)) x 0*fliplr(x(2:N+1))],N+1,1))).*Trap, 2))'/C_W), t_2_ref, s_2).*(1-x);
        options = optimoptions('fsolve','Display','none','OptimalityTolerance',...
            error_tolerance,'StepTolerance',error_tolerance);
        [G,fval,exitflag,output,jacobian] = fsolve(f,G,options); % use previous solution value as initial guess for nonlinear solve
        % store the relative error for analysis later
        error(jj,count) = delta*trapz(abs(fval));
        sol_norm(jj,count) = delta*trapz(G); % calculate L^1 norm of the solution
        %%
        figure(fig_count);
        %         % fancy plot can be commented out normally
        %         %area(0:delta:L,ones(1,N+1),'LineWidth',1.5,'FaceColor',[0 0.39 0]); % fancy plot
        %         %hold on; % fancy plot
        %         %area(0:delta:L,G,'LineWidth',1.5,'FaceColor',[0.565 0.933 0.565]); % fancy plot
        subplot(1,2,1), plot(0:delta:L,G,'-.b','LineWidth',2);
        xlim([0 L]); ylim([0 1]);
        xlabel('\Omega'); ylabel('Density');
        %         %legend('Forest','Grass');
        set(gca,'linewidth',1.25); set(gca,'FontSize',18);
        legend('IC','final soln');
        title(['\sigma = ',num2str(disp(count)),' and L^1 error = ',num2str(error(jj,count))]);
        z = eig(jacobian);
        principal_eig(jj,count) = max(real(z));
        subplot(1,2,2), scatter(real(z),imag(z));
        xlim([min(real(z))-0.05 max(real(z))+0.05]);
        ylim([min(imag(z))-0.025 max(imag(z))+0.025]);
        title(['Principal eigenvalue: ',num2str(max(real(z)))]);
        hold on;
        subplot(1,2,2), xline(0,'-.r','LineWidth',2);
    end
end
%% Plot the L^1 norm of the solution versus the dispersal parameter sigma
figure(100);
hold on;
L1_tolerance = 1e-4;
for i = 1:branches
    for j = 1:length(disp)
        if principal_eig(i,j) < 0 && error(i,j) < L1_tolerance
            scatter(disp(j),sol_norm(i,j),'b','filled');
        elseif principal_eig(i,j) < 0 && error(i,j) > L1_tolerance
            scatter(disp(j),sol_norm(i,j),'r','filled');
        elseif principal_eig(i,j) > 0 && error(i,j) > L1_tolerance
            scatter(disp(j),sol_norm(i,j),'r');
        else
            scatter(disp(j),sol_norm(i,j),'b');
        end
    end
end
hold on;
xlabel('\sigma');
ylabel('|G|');
xlim([0 0.35]);
ylim([0 1.05]);
yticks([0 0.2 0.4 0.6 0.8 1]);
set(gca,'linewidth',2);
set(gca,'FontSize',15);
%%
toc




