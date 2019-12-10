%% Plotting the Delta-Beta Ratio and the Sleep-Wake cycles for each Patient for each night

    patient = 'Patient_2';
    night = 'Night_2';
    cd(['/Users/bsevak/Documents/Merged Data_BF/Merged_Data/',patient,'/',night,'/']);

    load(['/Users/bsevak/Documents/Merged Data_BF/Merged_Data/',patient,'/',night,'/classifier_marker.mat']);

    figure;

    plot(Db);
    hold on;    
    plot(classifier_marker.wake_events,Db(classifier_marker.wake_events),'yo', 'MarkerSize', 5, 'MarkerFaceColor', 'y');
    hold on;
    plot(classifier_marker.sleep_events,Db(classifier_marker.sleep_events),'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'k')

    title(['Delta and Beta Ratio in Sleep wake states for ',patient,' ',night]);
    legend('Delta-Beta Ratio','Wake','Sleep');
    xlabel('Epoch Number (30 secs)');
    
    saveas(gcf,[patient,'_',night,'.fig']);
    
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [3 6 20 5];
    print([patient,'_',night],'-dpng','-r0')
    
   
 