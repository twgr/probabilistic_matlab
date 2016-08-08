function plot_lss_error(summary_lss,ground_truth,h1,h2,colour)
% Plots the errors for the LGSS model

squared_diffs = summary_squared_diffs(summary_lss,ground_truth);

plot_fields = {'Mean','Std_dev','Skewness','Excess_kurtosis'};

figure(h1);

for n=1:4;
    subplot(2,2,n);
    loglog(squared_diffs.conv_total.x.(plot_fields{n}),colour); 
    hold on;
    title(regexprep(plot_fields{n},'_',' '));
    xlabel('MCMC Iteration'); ylabel('Mean Squared Error');
    set(gca,'FontSize',32);
end

figure(h2);

for n=1:4;
    subplot(2,2,n);
    semilogy(mean(reshape(squared_diffs.pos_final.x.(plot_fields{n}),3,[]),1),colour); 
    hold on;
    title(regexprep(plot_fields{n},'_',' '));
    xlabel('State space time step t'); ylabel('Mean Squared Error');
    set(gca,'FontSize',32);
end