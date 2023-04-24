
close all
clear all


%% IF data directory

dir = '/Users/tim/Desktop/scripts_IF/example_process/output_2_IF_data/';

%% create a 96-well plate

Row=[1:8];
Col=[1:12];

dapiMatrix = cell(max(Row),max(Col));
nucfitcMatrix = cell(max(Row),max(Col));
nucCy5Matrix = cell(max(Row),max(Col));


%% Specify rows and cols for comparisons and load data

row=[1];  %[1:8] for 8 rows
col=[1];  %[1:12] for 12 columns

for n = row
    for i = col
        
        filename = [dir, num2str(n),'_', num2str(i),'_data'];
        load(filename);
        
        dapiMatrix{n,i} = intintdapi_allwells;
        nucfitcMatrix{n,i} = avgnucfitc_allwells;
        nucCy5Matrix{n,i} = avgnucCy5_allwells;
        
    end
end


%% DENSITY SCATTER PLOT

figure(1)
set(gcf, 'Renderer', 'painters');

dscatter([log(dapiMatrix{1,1})]', ([log(nucCy5Matrix{1,1})]'))
 xlim([12.7 14.5])
 ylim([5.7 8.6])
yline(6.85, '-r')
xline(13.65, '-r')
xlabel('DNA Content')
ylabel('Cyclin A2 (a.u.)')
title('MCF10A, 4h DMSO')

