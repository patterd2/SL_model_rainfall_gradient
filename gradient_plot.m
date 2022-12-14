P_fun = @(x, p_0, p_1) p_0 + p_1.*x/L;
P_fun_piece = @(x, p_0, p_1, slope,L) max(min( (p_0+p_1/2) + slope*(x-L/2)/L, p_0+p_1),p_0);
L = 1; % working on [0,L]
N = 300; % N+1 grid points
delta = L/N;
X = 0:delta:L;
%%
figure;
slope_vec = [1 2 3.5];
for i=1:length(slope_vec)
    plot(X,P_fun_piece(X,p_0,p_1,slope_vec(i),L),'-.','Linewidth',5);
    hold on;
    xlim([0 1]);
    ylim([0 1]);
    legend('slope = 1','slope = 2','slope = 3.5','FontSize',25);
end
xticks = 0:0.25:1;
set(gca, 'xtick', xticks);
yticks = 0:0.25:1;
set(gca, 'ytick', yticks);
ylabel('Precipitation Level');
grid on;
set(gca,'linewidth',2.5);
set(gca,'FontSize',25);