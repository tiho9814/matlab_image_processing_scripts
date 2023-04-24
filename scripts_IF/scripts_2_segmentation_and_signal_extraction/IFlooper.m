
%%  list all of your rows and columns and call the IF script


clear all

dir='/Users/tim/Desktop/scripts_IF/example_process/output_1_TIFF_files/';
mkdir(['/Users/tim/Desktop/scripts_IF/example_process/output_2_IF_data/'])

for row=[1]       %[1:8] for 8 rows
    for col=[1]   %[1:12] for 12 columns
        try
            
        [intintdapi_allwells, avgnucfitc_allwells, avgnucCy5_allwells, wellName]=IFfunction_TH_allCyto_allNuc(row, col, dir);
        
        save([ '/Users/tim/Desktop/scripts_IF/example_process/output_2_IF_data/' wellName, '_data'],'intintdapi_allwells', 'avgnucfitc_allwells', 'avgnucCy5_allwells')
        
        end
    end
end
