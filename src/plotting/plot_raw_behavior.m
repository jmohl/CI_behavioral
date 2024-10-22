%% raw behavior demo plots
%
% -------------------
% Jeff Mohl
% 7/24/19
% -------------------
%
% Description: generates plots of raw behavioral data for a few demo
% conditions. The purpose of these plots is to show that people and monkeys
% basically perform the task. There are 3 plots, which can be repeated for
% each subject or across subjects. These are:
% 1: unisensory localization accuracy for each target
% 2: saccade traces, with each frame plotted to show speed and fixation etc
% 3: saccade histograms for the same condition, maybe one each for
% integrated/segregated
function plot_raw_behavior(figpath,savefiles)
global aud_color vis_color

data_j = load('data\Juno_combined.mat');
data_j = data_j.tidy_data;

data_y = load('data\Yoko_combined.mat');
data_y = data_y.tidy_data;

data_j1 = load('data\Juno_AVD2_2018_04_11_tidy.mat'); %3_19_2, 4_11
data_j1 = data_j1.tidy_data;

data_y1 = load('data\Yoko_AVD2_2019_04_17_tidy.mat'); %4_17, 4_25 good examples but have non-0 fix
data_y1 = data_y1.tidy_data;

 data_H1 = load('data\H08_AVD2_2018_08_10_tidy.mat');
 data_H1 = data_H1.tidy_data;
% data_H1 = load('data\H07_AVD2_2018_08_09_tidy.mat');
% data_H1 = data_H1.tidy_data;

%combine human datasets
files_H = dir('data\*H*tidy*');
data_H = {};
for ind = 1:length(files_H)
    fname = files_H(ind).name;
    this_data = load(sprintf('data\\%s',fname));
    data_H = vertcat(data_H,this_data.tidy_data);
end

%rotate through subjects for each plots
 subjects = {data_j,data_y,data_H,data_j1,data_y1,data_H1};
 subject_id = {'Juno','Yoko','H','Juno1day','Yoko1day','H07'};
%   subjects = {data_j1,data_y1};
%   subject_id = {'Juno1day','Yoko1day'};
% subjects = {data_H1};
% subject_id = {'H07'};
horz_only = 0; %option to plot horizontal component of eye trace over time

%% run all plots on all subjects
for ind = 1:length(subjects)
    figure
    this_data = subjects{ind};
    this_data = this_data(logical(this_data.valid_tr),:); %only include valid trials
    
    %% plot 1, unisensory localization error bar plots
    %bars represent STD for that day, to show that variance is higher on aud
    %trials
    jitter = .25;
    A_data = this_data(strcmp(this_data.trial_type,'A'),:);
    V_data = this_data(strcmp(this_data.trial_type,'V') & abs(this_data.V_tar) ~= 30,:);
    
    [gA,glabA] = findgroups(A_data.A_tar);
    [gV,glabV] = findgroups(V_data.V_tar);
    A_sac = vertcat(A_data.valid_endpoints{:});
    A_sac = A_sac(:,1);
    V_sac = vertcat(V_data.valid_endpoints{:});
    V_sac = V_sac(:,1);
    
    mean_A = splitapply(@mean,A_sac,gA);
    mean_V = splitapply(@mean,V_sac,gV);
    std_A = splitapply(@std,A_sac,gA);
    std_V = splitapply(@std,V_sac,gV);
    
    subplot(1,3,1)
    errorbar(glabA+jitter,mean_A,std_A,'Color',aud_color,'LineWidth',1)
    grid on
    hold on
    errorbar(glabV-jitter,mean_V,std_V,'Color',vis_color,'LineWidth',1)
    labels = [-24 -18 -12 -6 0 6 12 18 24];
    xticks(labels)
    xticklabels(labels)
    yticks(labels)
    yticklabels(labels)
    xlim([-30 30])
    ylim([-30 30])
    ylabel('Saccade endpoint')
    xlabel('Target Location')
    legend('Auditory trials','Visual trials','Location','Best')
    title(sprintf('Accuracy of single target saccades: %s',subject_id{ind}),'Interpreter','none')
    
    %% plot 2, saccade eye traces for an example AV condition with two saccades
    if strfind(subject_id{ind},'H')
        ex_tars = [12, -12];
    else
        ex_tars = [6, -12];
    end
    AV_ex_data = this_data(this_data.A_tar == ex_tars(1) & this_data.V_tar == ex_tars(2),:);
    subplot(1,3,2)
    %if the number of trials is really large, need to subsample down to a
    %manageable number
    if height(AV_ex_data) > 30
        rng('default')
        AV_ex_data = AV_ex_data(randsample(height(AV_ex_data),30),:);
    end
    
    if horz_only
        hold on
        ref_tar_A = plot([-300,1700],[ex_tars(1), ex_tars(1)],'--r');
        ref_tar_V = plot([-300,1700],[ex_tars(2), ex_tars(2)],'--b');
        ref_go = plot([0 0],[-20,20],'k-','LineWidth',.5);
        for tr = 1:height(AV_ex_data)
            
            eyedata = AV_ex_data(tr,:).eyedata{:};
            go_time = AV_ex_data(tr,:).go_time;
            end_time = AV_ex_data(tr,:).end_time;
            start_int = go_time - 300;
            time_vector = -300:2:(end_time-go_time);
            plot(time_vector,eyedata(start_int/2:end_time/2,1),'Color',[.75 .75 .75],'LineWidth',.2)
            
        end
        xlim([-300,1700])
        ylim([-20,20])
        title(sprintf('%s, %dA %dV',subject_id{ind},ex_tars))
        xlabel('time(ms)')
        ylabel('Vertical eye position (deg)')
        legend([ref_tar_A,ref_tar_V, ref_go],'A target','V target','Go cue')
        
    else
        
        hold on
        ref_tar_A = plot([ex_tars(1), ex_tars(1)],[-30,30],'--','color',aud_color);
        ref_tar_V = plot([ex_tars(2), ex_tars(2)],[-30,30],'--', 'color',vis_color);
        
        for tr = 1:height(AV_ex_data)
            %subsample eye data, 1 dot per x*2 ms, remember sampling rate = 500hz
            subsamp = 2;
            eyedata = AV_ex_data(tr,:).eyedata{:};
            %only get eye data between go cue and end
            eyedata = eyedata(AV_ex_data(tr,:).go_time/2:AV_ex_data(tr,:).end_time/2,:); %divide by 2 to make in eye time
            eyedata = eyedata(1:subsamp:end,:);
            %             plot(eyedata(:,1),eyedata(:,2),'Color',[.75 .75 .75],'LineWidth',.2)
            plot(eyedata(:,1),eyedata(:,2),'k.','MarkerSize',.05)
        end
        title('Double target saccades')
        xlabel('Horizontal eye position (deg)')
        ylabel('Vertical eye position (deg)')
        legend([ref_tar_A,ref_tar_V],'A target','V target' )
        
        %plot visual and aud ref lines
        xlim([-30,30])
        ylim([-30,30])
    end
    
    %% plot3, histogram of saccade endpoints.
    % This is essentially copying plot_localization but I am working with the
    % data in a different format. Probably worth considering just using
    % plot_localization instead so that it exactly matches

    AV_sac_A = vertcat(AV_ex_data.A_endpoints{:});
    AV_sac_V = vertcat(AV_ex_data.V_endpoints{:});

       
    subplot(1,3,3)
    hold on
    ref_tar_A = plot([ex_tars(1), ex_tars(1)],[0,.35],'--','color',aud_color,'LineWidth',2);
    ref_tar_V = plot([ex_tars(2), ex_tars(2)],[0,.35],'--','color',vis_color,'LineWidth',2);
  
    A_hist = histogram(AV_sac_A(:,1),-30:30,'Normalization','Probability','FaceColor',aud_color,'FaceAlpha',1);
    V_hist = histogram(AV_sac_V(:,1),-30:30,'Normalization','Probability','FaceColor',vis_color,'FaceAlpha',1);
    title('Saccade endpoints')
    xlabel('Horizontal eye position (deg)')
    ylabel('p in bin')
    legend([A_hist,V_hist,ref_tar_A,ref_tar_V],'Labeled Auditory', 'Labeled Visual','A target','V target')

    % below, I toyed with showing the unimodal distribution from single
    % target trials and/or means taken from these distributions, but
    % decided it was too cluttered.
%plot unimodal dist
%     A_ex_data = A_data(A_data.A_tar == ex_tars(1),:);
%     V_ex_data = V_data(V_data.V_tar == ex_tars(2),:);

%     A_sac = vertcat(A_ex_data.valid_endpoints{:});
%     V_sac = vertcat(V_ex_data.valid_endpoints{:});
%     xrange = -50:50;
%     A_sac_norm = histcounts(A_sac(:,1),xrange)/length(A_sac(:,1));
%     V_sac_norm = histcounts(V_sac(:,1),xrange)/length(V_sac(:,1));
%     mean_A = mean(A_sac(:,1));
%     mean_V = mean(V_sac(:,1));

%     plot(-49.5:49.5,A_sac_norm,'--','Color',[1 0 0 .5]);
%     plot(-49.5:49.5,V_sac_norm,'--','Color',[0 0 1 .5]);
  % plot mean indicator - triangle
  %     AV_mean_A = mean(AV_sac_A(:,1));
%     AV_mean_V = mean(AV_sac_V(:,1));
%     max_p = .25;
%     plot(mean_A,max_p,'v','MarkerSize',8, 'MarkerEdgeColor',aud_color)
%     plot(mean_V,max_p,'v','MarkerSize',8, 'MarkerEdgeColor',vis_color)
%     plot(AV_mean_A,max_p,'v','MarkerSize',8,'MarkerEdgeColor',aud_color,'MarkerFaceColor',aud_color)
%     plot(AV_mean_V,max_p,'v','MarkerSize',8,'MarkerEdgeColor',vis_color,'MarkerFaceColor',vis_color)
%     
%     % do ttest between the two distributions
%     [HA,pA] = ttest(AV_sac_A(:,1)-mean_A);
%     [HV,pV] = ttest(AV_sac_V(:,1)-mean_V);
%     if HA
%         text(mean_A,max_p+.02,sprintf('%d',pA))
%     end
%     if HV
%         text(mean_V,max_p+.02,sprintf('%d',pV))
%     end    
%     

    %% save out figure
    set(gcf,'Position',[25,50,1400,350])
    if savefiles
        saveas(gcf,sprintf('%s\\%s_rawdata',figpath,subject_id{ind}),'svg');
        saveas(gcf,sprintf('%s\\%s_rawdata',figpath,subject_id{ind}),'png');
    end
    

end





