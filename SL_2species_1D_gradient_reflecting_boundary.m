%% 1D Simulation of spatial savanna model with precipitation
% Explicit Euler in time, 2D trapezoidal rule in space
% the boundaries are reflecting
% NB fire is modelled as the integral of phi of G
tic
%% Set options for plots/movies
close all
fprintf('\n');
DO_MOVIE = 0; % i.e. write movie to avi file for playback, else just plots
SPACE_TIME_PLOT = 0;
FINAL_PLOT = 0;
%% Numerical method parameters
L = 1; % working on [0,L]
N = 500; % N+1 grid points
delta = L/N;  % spatial discretization parameter
h = 0.05; % time discretisation parameter
n = 10000; % number of time steps
tau = (n-1)*h; % simulations time domain is [0,tau]
%% Function definitions
P_fun = @(x, p_0, p_1) p_0 + p_1.*x;
alpha_p = @(alpha, alpha_s, p) alpha + alpha_s.*p;
J_F_fun = @(x, a, sigma_F) exp( -( (a-x).^2)/(2*sigma_F^2) )/(sqrt(2*pi*sigma_F^2)); % 1./(pi*sigma_F*(1+((x-a).^2)/sigma_F^2));
W_fun = @(x, a, sigma_W) exp( -( (a-x).^2 )/(2*sigma_W^2) )/(sqrt(2*pi*sigma_F^2)); % 1./(pi*sigma_W*(1+((x-a).^2)/sigma_W^2));
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
s_2 = 0.05;% s_2's standard value is 0.05

%linestyles = ['-',':','-.','-',':','-.'];
fill_color = ['r','k','b','c','m','g'];
disp = 0.01;
ICs = 2;
sol_norm = zeros(length(disp),length(ICs));
relative_error = zeros(length(disp),n-1);

for IC = ICs % 1 for high grass IC, 2 for low grass IC, 4 for mixed spatial stripes between the two species, 
    % 3 for a sigmoidal transition, 6+ for random IC
    for count = 1:length(disp)
        
        sigma_F = disp(count); % seed dispersal radius forest trees
        sigma_W = disp(count);%disp(count); % fire spread radius
        
        %% Set up the initial distributions of the cover types on the grid
        % each row is one time step of the simulation
        % solution is a block matrix [LB; SOL; RB] where LB and RB are fixed
        % boundary conditions corresponding to a "reflecting-type boundary"
        if IC == 1
            G0 = 1-0.05*(1+cos(12*pi*(0:delta:L)));
        elseif IC == 2
            G0 = 0.05*(1+cos(12*pi*(0:delta:L)));
        elseif IC == 4
            G0 = 0.5*(1+cos(12*pi*(0:delta:L)));
        elseif IC == 3
            G0 = phi(0.95, 0.01, 0:delta:L, 0.5, 0.05);
        elseif IC == 5
            G0 = 0.5*(1+cos(24*pi*(0:delta:L)));
        else
            %G0 = rand(1,N+1);
        end
        LB_G = fliplr(G0(1,2:end));
        RB_G = fliplr(G0(1,1:end-1));
        G = [LB_G G0 RB_G];
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
            integrand = tempW(i,:).*G(1,:);
            E(1,i) = sum(integrand.*Trap)*delta;
        end
        C_W = max(temp_normalise);
        E(1,:) = E(1,:)/C_W;
        %% Compute the birth and mortality matrices as a function of rainfall
        P_grad = P_fun(X,p_0,p_1); % compute the rainfall gradient along the x-axis
        alpha_grad = alpha_p(alpha, alpha_s, P_grad);
        %% preallocaate some temp variables for efficiency
        temp1 = ones(1,N+1);
        tempF = ones(N+1,3*N+1);
        temp_normalise_F = ones(1,N+1);
        %% pre-calculate the 4D convolution matrices
        for k = 1:N+1
            tempF(k,:) = J_F_fun([X_L X X_R],(k-1)*delta,sigma_F);
            temp_normalise_F(1,k) = sum(tempF(k,:))*delta;
        end
        C_F = max(temp_normalise_F);
        %% The numerical scheme
        for i = 2:n
            % compute convolutions for this time step
            progressbar(i,n);
            for k = 1:N+1
                integrand1 = tempF(k,:).*(1-G(i-1,:));
                temp1(1,k) = sum(integrand1.*Trap)*delta;
            end
            temp1 = temp1/(C_F);
            G(i,(N+1):2*N+1) = G(i-1,(N+1):2*N+1) + h*( - alpha_grad.*temp1.*G(i-1,(N+1):2*N+1) + ...
                phi(f_0_ref, f_1_ref, E(i-1,:), t_2_ref, s_2).*(1-G(i-1,(N+1):2*N+1)) );
            % now need to update the extended parts of the solution
            G(i,1:N) = fliplr(G(i,N+2:2*N+1));
            G(i,2*N+2:3*N+1) = fliplr(G(i,(N+1):2*N));
            
            for k = 1:N+1
                integrand = tempW(k,:).*( G(i,:) );
                E(i,k) = delta*sum(integrand.*Trap)./C_W;
            end
            % For numerical stability take max with zero and min with one
            % add error message here later
            G(i,:) = min(max(G(i,:),0),1);
            E(i,:) = min(max(E(i,:),0),1);
            % store the relative error for analysis later
            relative_error(count,i-1) = max(abs(G(i,:)-G(i-1,:)))/max(abs(G(i-1,:)));
        end
        sol_norm(count,IC) = delta*sum(G(end,(N+1):2*N+1)); % calculate L^1 norm of the solution
        %fprintf('\n');
        %fprintf('Relative error in the discrete sup norm is:\n');
        %fprintf(['G: ',num2str(max(abs(G(n,:)-G(n-1,:)))/max(abs(G(n-1,:)))),'\n']);
        %%
        if FINAL_PLOT
            figure(1);
            %f.Name = ['Simulation time: t = ', num2str((j-1)*h)];
            %subplot(2,2,IC), plot(X,G(end,N+1:2*N+1),'r','LineWidth',2,'LineStyle',linestyles(count));
            subplot(2,2,1), area(0:delta:L,ones(1,N+1),'LineWidth',1.5,'FaceColor',[0 0.39 0]);
            hold on;
            subplot(2,2,1), area(0:delta:L,G(20,N+1:2*N+1),'LineWidth',1.5,'FaceColor',[0.565 0.933 0.565]);
            title('t = 1');
            xlim([0 L]);
            ylim([0 1]);
            %
            subplot(2,2,2), area(0:delta:L,ones(1,N+1),'LineWidth',1.5,'FaceColor',[0 0.39 0]);
            hold on;
            subplot(2,2,2), area(0:delta:L,G(100,N+1:2*N+1),'LineWidth',1.5,'FaceColor',[0.565 0.933 0.565]);
            title('t = 5');
            xlim([0 L]);
            ylim([0 1]);
            %
            subplot(2,2,3), area(0:delta:L,ones(1,N+1),'LineWidth',1.5,'FaceColor',[0 0.39 0]);
            hold on;
            subplot(2,2,3), area(0:delta:L,G(200,N+1:2*N+1),'LineWidth',1.5,'FaceColor',[0.565 0.933 0.565]);
            title('t = 10');
            xlim([0 L]);
            ylim([0 1]);
            %
            subplot(2,2,4), area(0:delta:L,ones(1,N+1),'LineWidth',1.5,'FaceColor',[0 0.39 0]);
            hold on;
            subplot(2,2,4), area(0:delta:L,G(end,N+1:2*N+1),'LineWidth',1.5,'FaceColor',[0.565 0.933 0.565]);
            title('Final time');
            xlim([0 L]);
            ylim([0 1]);
            legend('Forest','Grass','Location','SouthWest');
            set(gca,'linewidth',1.25);
            set(gca,'FontSize',18);
            %             if IC == 1
            %                 dim = [.24 .615 .3 .3];
            %                 str = sprintf('\\bf\\sigma = %.2f', disp(count));
            %                 annotation('textbox',dim,'String',str,'FitBoxToText','on','LineWidth',1.25,'FontSize',11,'Margin',2.5,'HorizontalAlignment','center');
            %             elseif IC == 2
            %                 dim = [.68 .615 .3 .3];
            %                 str = sprintf('\\bf\\sigma = %.2f', disp(count));
            %                 annotation('textbox',dim,'String',str,'FitBoxToText','on','LineWidth',1.25,'FontSize',11,'Margin',2.5,'HorizontalAlignment','center');
            %             elseif IC == 3
            %                 dim = [.68 .135 .3 .3];
            %                 str = sprintf('\\bf\\sigma = %.2f', disp(count));
            %                 annotation('textbox',dim,'String',str,'FitBoxToText','on','LineWidth',1.25,'FontSize',11,'Margin',2.5,'HorizontalAlignment','center');
            %             elseif IC == 4
            %                 dim = [.68 .135 .3 .3];
            %                 str = sprintf('\\bf\\sigma = %.2f', disp(count));
            %                 annotation('textbox',dim,'String',str,'FitBoxToText','on','LineWidth',1.25,'FontSize',11,'Margin',2.5,'HorizontalAlignment','center');
            %             end
            xlabel('\Omega');
            ylabel('Density');
            set(gca,'linewidth',1.25);
            set(gca,'FontSize',12);
            
            hold on;
        end
        
    end
end
%% Space time plot of the solution
fST = figure;
fST.Name = 'Evolution over time';
%subplot(2,2,1),
h1 = pcolor(G(:,N+1:2*N+1));
% custom_map = [
%     linspace(1,0,100)' linspace(1,0.5,100)' linspace(1,0,100)'];
% colormap(custom_map);
shading interp
%title('Grass');
c = colorbar;
c.LineWidth = 1.5;
set(h1, 'EdgeColor', 'none');
ylabel('Time');
caxis([0,1])
xticks([1 floor((N+1)/2) floor((N+1))]);
xticklabels({num2str(0), num2str(L/2), num2str(L)});
yticks([1 floor(n/2) floor(n)]);
yticklabels({num2str(0), num2str(n*h/2), num2str(n*h)});
xlabel('\Omega');
ylabel('Time');
set(gca,'linewidth',4);
set(gca,'FontSize',20);

%% Plot the L^1 norm of the solution versus the dispersal parameter
% figure(5);
% for i = 1:3
%     hold on;
%     scatter(disp,sol_norm(:,i),'filled',fill_color(i));
%     hold on;
% end
% xlabel('\sigma');
% ylabel('|G|');
% xlim([0.1 0.3])
% ylim([0 1])
% yticks([0 0.2 0.4 0.6 0.8 1]);
% set(gca,'linewidth',2);
% set(gca,'FontSize',15);

%%
toc
% output the relative errors
relative_error(:,end)

%% stacked plots 


