
function umich_ott

umich_ott=[49 46 51 42 45 43 50 51 56 46 50 64 44 47 42 62 45 61 56 56 41 43 ...
43 43 57 42 41 43 45 49 59 62 48 42 48 49 43 51 44 44 63 63 46 59 ...
47 50 42 42 42 63 63 44 71 50 50 47 55 69 70 59 46 57 74 58 46 61 ...
50 58 49 62 45 44 66 49 48 56 57 59 50 45 70 62 49 55 71 49 45 45 69];

boxplot(umich_ott);
xlabel('UMich','FontSize',14);
ylabel('OTT (days)','FontSize',14);
text(1.1,50,['Median OTT: ',num2str(median(umich_ott))],...
    'FontSize',14);

end