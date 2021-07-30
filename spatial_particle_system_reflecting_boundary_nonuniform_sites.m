%% Matrix simulation of Poisson Forest-Grass Markov Chain
clear all;

% RNG seed
%rng(1546793)

Times=[];
DO_PLOT=0;
No_Abs=1;

%% Spatial Parameters
L = 1; % working on [0,L]
% working with a uniform site distribution

%% Model Parameters
alpha=0.5; 
alpha_s = 1.25; % slope for the gradient, [alpha, alpha + alpha_s] is the range

t2=0.4;
f0=0.1;
f1=0.9;
s2=0.05;

%% dispersal and fire kernels
J_F_fun = @(x, a, sigma_F) 1./(pi*sigma_F*(1+((x-a).^2)/sigma_F^2));%exp( -( (a-x).^2)/(2*sigma_F^2) )/sqrt(2*pi*sigma_F^2);
W_fun = @(x, a, sigma_W) 1./(pi*sigma_W*(1+((x-a).^2)/sigma_W^2));%exp( -( (a-x).^2 )/(2*sigma_W^2) )/sqrt(2*pi*sigma_W^2);
% precipitation gradient impact on growth
alpha_p = @(alpha, alpha_s, p) alpha + alpha_s.*p;
phi=@(x) (f0+(f1-f0)./(1+exp(-(x-t2)/s2))); % fire sigmoid
phi2 = @(f_0_fun, f_1_fun, g, t_2_fun, s_2_fun) f_0_fun + (f_1_fun - f_0_fun)./(1 + exp(-(g-t_2_fun)/s_2_fun));

sigma_F = 0.025; % seed dispersal radius forest trees
sigma_W = sigma_F; % fire spread radius

sites = 1000;

P=1;                % Number of Patches
N=sites*ones(1,P);    % Number of Sites / Patch
NTot=sum(N);        % Total number of sites

T=500; % max time
dt=0.01; % for interpolating uniformly in time later
t0=0;

MC=1;

J=alpha;
W=1;
Ntimes=length(t0:dt:(T-dt));

tic

nbins=100;
Histo=zeros(nbins,MC);
p_grass = 0.9;
p_no_grass = 1 - p_grass;
%% This block of code runs one simulation via the Gillespie algorithm
for iteration=1:MC
    Solution=zeros(NTot,1);
    Locations=sort(rand(sites,1)); % inverse transform method for random site locations
    alpha_grad = alpha_p(alpha, alpha_s, Locations);
    % histogram(Locations, 50); % check the dsitribution of the sites
    [A,B]=meshgrid(Locations);
    [A_left,B_left]=meshgrid(Locations-L);
    [A_right,B_right]=meshgrid(Locations+L);
    
    % need to calculate the interactions from [-L,0] and [L,2L] to the
    % solution on [0,L]
    J_Mat = (J_F_fun(A,B,sigma_F) + fliplr(J_F_fun(A,B_left,sigma_F)) + fliplr(J_F_fun(A,B_right,sigma_F)))/sites;
    W_Mat = (W_fun(A,B,sigma_W) + fliplr(W_fun(A,B_left,sigma_W)) + fliplr(W_fun(A,B_right,sigma_W)))/sites;
    
    % check that the normalization of the kernels is correct
    %[X,Y] = meshgrid(0:1/sites:L);
    %J_Mat_regular = J_F_fun(X,Y,sigma_F)/sites;
    %plot(max(sum(J_Mat_regular)))
    
    % set the initial conditions
    Solution(:,1) = rand(NTot,1) < p_grass; %phi2(0.95, 0.05, Locations, 0.5, 0.05);
    Times = [0];
    
    k=0;
    t=0;
    while (t<T)
        k=k+1;
        BirthRates=alpha_grad.*((J_Mat*(1-Solution(:,k)))).*Solution(:,k);
        DeathRates=phi(W_Mat*Solution(:,k)).*(1-Solution(:,k));
        %BirthRates=alpha*((J_Mat*(1-Solution(:,k)))).*Solution(:,k);
        %DeathRates=phi(W_Mat*Solution(:,k)).*(1-Solution(:,k));
        
        totalIntensity=sum(BirthRates+DeathRates);
        
        NextEvent=-log(1-rand())/totalIntensity;
        
        t=t+NextEvent;
        Times(end+1)=t;
        
        CDF=cumsum(BirthRates+DeathRates)/totalIntensity;
        U=rand();
        i=1;
        while U>CDF(i)
            i=i+1;
        end
        Solution(:,k+1)=Solution(:,k);
        Solution(i,k+1)=1-Solution(i,k);
    end
    k=k+1;
    Times(end)=T;
    %Solution(:,k+1)=Solution(:,k);
end
toc

%% Create evenly spaced solution data for plots

% linearly interpolate solution values to put simulation on true timescale
sol_new = ones(sites, length(0:dt:Times(end)));
for i = 1:sites
    sol_new(i,:) = interp1(Times,Solution(i,:),0:dt:Times(end));
end

% need to put the spatial dimension on a linear grid as well
dx = 0.005; % spatial mesh
sol_new2 = ones(length(0:dx:Locations(end)),length(0:dt:Times(end)));
for i=1:length(0:dt:Times(end))
    sol_new2(:,i) = interp1(Locations,sol_new(:,i),0:dx:Locations(end));
end

%%
% Space-time plot
figure(4);
imagesc(flipud(sol_new2'));
custom_map = [1 1 1
    0 0.5 0];
colormap(custom_map);
caxis([0 1]);
xlabel('\Omega');
ylabel('Time');
dx_max = 1/dx;
dt_max = T/dt;
xticks([1 floor((dx_max+1)/4) floor((dx_max+1)/2) floor(3*(dx_max+1)/4) floor((dx_max))]);
xticklabels({num2str(0), num2str(L/4), num2str(L/2), num2str(3*L/4), num2str(L)});
yticks([1 floor(dt_max/4) floor(dt_max/2) floor(3*dt_max/4) floor(dt_max)]);
yticklabels({num2str(T), num2str(T*3/4), num2str(T/2), num2str(T/4), num2str(0)});
set(gca,'FontSize',22);
set(gca,'linewidth',1.75);

% plot the average proportion of grass in the system 
% should be able to roughly observe the Maxwell point

%%
% Plot of cross section of the final solution
figure(2);
subplot(2,2,1), plot(mean(sol_new2(:,floor(0.2*length(0:dt:Times(end))):floor(0.25*length(0:dt:Times(end)))),2),'LineWidth',2);
xticks([1 floor((dx_max+1)/4) floor((dx_max+1)/2) floor(3*(dx_max+1)/4) floor((dx_max))]);
xticklabels({num2str(0), num2str(L/4), num2str(L/2), num2str(3*L/4), num2str(L)});
ylim([0 1]);
title(sprintf('T = %.0f', 0.25*T));
xlabel('\Omega');
ylabel('Grass density');
set(gca,'linewidth',1.25);
set(gca,'FontSize',12);


subplot(2,2,2), plot(mean(sol_new2(:,floor(0.45*length(0:dt:Times(end))):floor(0.5*length(0:dt:Times(end)))),2),'LineWidth',2);
xticks([1 floor((dx_max+1)/4) floor((dx_max+1)/2) floor(3*(dx_max+1)/4) floor((dx_max))]);
xticklabels({num2str(0), num2str(L/4), num2str(L/2), num2str(3*L/4), num2str(L)});
ylim([0 1]);
title(sprintf('T = %.0f', 0.5*T));
xlabel('\Omega');
set(gca,'linewidth',1.25);
set(gca,'FontSize',12);


subplot(2,2,3), plot(mean(sol_new2(:,floor(0.7*length(0:dt:Times(end))):floor(0.75*length(0:dt:Times(end)))),2),'LineWidth',2);
xticks([1 floor((dx_max+1)/4) floor((dx_max+1)/2) floor(3*(dx_max+1)/4) floor((dx_max))]);
xticklabels({num2str(0), num2str(L/4), num2str(L/2), num2str(3*L/4), num2str(L)});
ylim([0 1]);
title(sprintf('T = %.0f', 0.75*T));
xlabel('\Omega');
ylabel('Grass density');
set(gca,'linewidth',1.25);
set(gca,'FontSize',12);


subplot(2,2,4), plot(mean(sol_new2(:,floor(0.95*length(0:dt:Times(end))):end),2),'LineWidth',2);
xticks([1 floor((dx_max+1)/4) floor((dx_max+1)/2) floor(3*(dx_max+1)/4) floor((dx_max))]);
xticklabels({num2str(0), num2str(L/4), num2str(L/2), num2str(3*L/4), num2str(L)});
ylim([0 1]);
title(sprintf('T = %.0f', T));
xlabel('\Omega');
set(gca,'linewidth',1.25);
set(gca,'FontSize',12);

%plot(mean(sol_new2(:,floor(0.975*length(0:dt:Times(end))):end),2),'LineWidth',2);
%scatter(Locations,Solution(:,end));


%%
% % Plot the interpolated (in time) solution as a 3D surface
% figure(3)
% sol_smooth = zeros(length(0:dx:Locations(end)), length(0:dt:Times(end)));
% for i = 1:length(Times)
%    sol_smooth(:,i) = mean(sol_new2(:,floor(0.9*length(0:dt:Times(end))):end),2); 
% end
% [X,Y] = meshgrid(0:dt:T,dx:dx:1);
% surf(X,Y,sol_smooth);
% shading interp
% colormap(custom_map);


