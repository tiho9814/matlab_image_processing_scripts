clear all;

% folder containing necessary functions for transforming tracedata
% function folder consists of files originally created by Mingwei Min
addpath ('Mingwei_timeslicer');

%% Common parameters to set
% Starting values correspond to the example tracedata file titled 'tracedata_1_1_1.mat'
% Example tracedata consists of DMSO-treated MCF10A cells expressing DHB-mVenus
% Check tracedata to determine value for each parameter

nuc_col     = 7;        % col in tracedata for cdk2 nuclear intensity
cyto_col    = 9;        % col in tracedata for cdk2 cyto intensity
other_col   = 5;        % protein of interest - for instance H2B/mCherry here
framerate   = 5;        % frames per hour
datadir     = '/Users/tim/Desktop/scripts_tracedata_plotting/';    % enter directory containing tracedata here
savedir     = '/Users/tim/Desktop/scripts_tracedata_plotting/';    % enter directory for saving cellinfo files here
movie_leng  = 260;      % total number of frames in tracedata
sitemat     = 1;      % total number of sites per well imaged
plot_time   = [-10 20]; % time range to plot
plot_signal = [30 180]; % signal range of the protein of interest to plot
drug_frame  = 127;      % frame in the tracedata where drug was added
windowmat = [4 10];    % window of POI2 relative to initial POI to slice for
max_trace   = 100;      % maximum number of traces to plot per condition
filter_col  = 9;
matr1=[1];            % columns to be analyzed together as a set


%% Loop through tracedata files of interest

tic

for kk=1  
    colmat=matr1(kk,:);
    for rowmat = [1]     %[1:8] for analyzing 8 different conditions laid out in rows
        if rowmat == 1
            drug = 'Treatment 1';
        elseif rowmat == 2
            drug = 'Treatment 2';
        elseif rowmat == 3
            drug = 'Treatment 3';
        elseif rowmat == 4
            drug = 'Treatment 4';
        elseif rowmat == 5
            drug = 'Treatment 5';
        elseif rowmat == 6
            drug = 'Treatment 6';
        elseif rowmat == 7
            drug = 'Treatment 7';
        elseif rowmat == 8
            drug = 'Treatment 8';
        end
        
        for aa=1
            grab_time   = windowmat(aa,:);

       %% grabbing cell traces of interest using MM functions
       
            [cdk2_trace, other_trace] = grab_align('colmat', colmat, ...
                'rowmat', rowmat, ...
                'sitemat', sitemat, ...
                'datapath', datadir, ...
                'cyto_col', cyto_col, ...
                'nuc_col', nuc_col, ...
                'other_col', other_col, ...
                'framerate', framerate, ...
                'movie_leng', movie_leng, ...
                'POI_name', 'mitosis', ...       % point of interest 1: traces are algined to this point. current options: mitosis, drug_addition (same as timeslicer), R_point, mitosis_skip
                'grab_cdk2state','cdk2inc', ...  % current options: all, cdk2inc, cdk2low, cdk2emerge
                'grab_time', grab_time,...
                'POI2_name', 'drug_addition',... % point of interest 2: with POI1 to creat a time interval of interest. Same options as POI1
                'drug_frame', drug_frame,...
                'filter_col', filter_col);
            

       %% Ploting single-cell and median CDK2 tracedata
            
         time   = (-movie_leng:movie_leng)/framerate;
         colors = [0 0 0;  %%% 000 is black, 010 is green, 001 is blue, 100 is red
                   0 1 0;
                   0 0 1;
                   1 0 0];
            
         % pause
         figure, hold on
            
         %% ploting single CDK2 traces
         subplot(2,2,1)
            
            patch([grab_time grab_time(2) grab_time(1)],[0.21 0.21 2 2],[0.85 0.85 0.85],'LineStyle','none');
            for i = 1:min(size(cdk2_trace,1),max_trace)
                line(time,cdk2_trace(i,:),'color','b');  % 'b' for blue
                hold on
            end

            xlim(plot_time);
            ylim([0.0 2.5]);
            xlabel('Time relative to anaphase (h)')
            ylabel('DHB C/N ratio')
            title(num2str(i)) % records number of traces plotted in title

         %% ploting median CDK2 trace that corresponds to single-cell plot
         subplot(2,2,2)

            patch([grab_time grab_time(2) grab_time(1)],[0.21 0.21 2 2],[0.85 0.85 0.85],'LineStyle','none');
            ciplot_pro2(cdk2_trace,time,colors(3,:));
            
            xlim(plot_time);
            ylim([0.2 2.0]);
            xlabel('Time relative to anaphase (h)')
            ylabel('Median CDK2 activity')
            set(gcf,'Color','w');

        end
    end
end

toc
