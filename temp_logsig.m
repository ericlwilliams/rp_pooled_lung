function temp_logsig

fig_loc = 'Z:/elw/MATLAB/meta_analy/slides/figures/latest/';

 n = -5:0.05:5;
a = logsig(n);
fig_sig=figure;
plot(n,a,'k','LineWidth',2)
ylabel('Probability of Complication','FontSize',14);

set(gca,'FontSize',12);
set(gca,'XTick',-5:1:5);
%set(gca,'XTickLabel',0:1:10);
set(gca,'XTickLabel',[]);
ylim([0 1]);

set(fig_sig,'Color','w');
export_fig(fig_sig,[fig_loc,'sig_response'],'-jpg');
disp(['Saving ',fig_loc,'sig_response.pdf...']);

end
