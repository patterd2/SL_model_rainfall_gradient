figure;

slope_vec = [1 2 3.5];
for i=1:length(slope_vec)
    plot(X,P_fun_piece(X,p_0,p_1,slope_vec(i)),'-.','Linewidth',3);
    hold on;
    xlim([0 1]);
    ylim([0 1]);
    legend('slope = 1','slope = 2','slope = 3.5','FontSize',22);
end
xticks = 0:0.25:1;
set(gca, 'xtick', xticks);
yticks = 0:0.25:1;
set(gca, 'ytick', yticks);
ylabel('Precipitation Level');
grid on;
set(gca,'linewidth',2);
set(gca,'FontSize',20);