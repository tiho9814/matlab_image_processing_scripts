
%% This code segments  nuclei in dapi, written for 3-color images

function[intintdapi_allwells, avgnucfitc_allwells, avgnucCy5_allwells, wellName]= IFfunction_TH_allCyto_allNuc(row, col, dir)

row 
col

intintdapi_allwells=[];
avgnucfitc_allwells=[];
avgnucCy5_allwells=[];  
 

for site=[1:6] %number of sites per  well
    
   wellName=[num2str(row),'_', num2str(col)];
    
   %load images
   
    dapi=single(imread([dir, num2str(row), '_', num2str(col), '_', num2str(site), '_BFP_1.tif']));  %the part in purple in quotes has to match file name
    fitc=single(imread([dir, num2str(row), '_', num2str(col), '_', num2str(site), '_FITC_1.tif']));
    Cy5=single(imread([dir, num2str(row), '_', num2str(col), '_', num2str(site), '_Cy5_1.tif']));
    
   dapimask=getdapimask(dapi, 11);       % the number after the comma is the nuclear RADIUS (not diameter); whole numbers only
   figure(1); imshowc(dapimask, dapi, 1) % plot nuclear mask to check quality of segmentation
   numberednuc=bwlabel(dapimask);        % give numbers to each nucleus  
   nucxypos=regionprops(numberednuc, 'PixelIdxList', 'Centroid', 'Area');
    
    for i=1:length(nucxypos)  %for each cell in the image
        
        sumnucdapi(i)=sum(dapi(nucxypos(i).PixelIdxList));
        avgnucfitc(i)=mean(fitc(nucxypos(i).PixelIdxList));
        avgnucCy5(i)=mean(Cy5(nucxypos(i).PixelIdxList));
                
       %get the integrated intensity of dapi stain to get DNA content, as opposed to the mean
        center(i,:)=nucxypos(i).Centroid;
    end
    
   %Put all wells together
    intintdapi_allwells=[intintdapi_allwells, sumnucdapi]; 
    avgnucfitc_allwells=[avgnucfitc_allwells, avgnucfitc];
    avgnucCy5_allwells=[avgnucCy5_allwells, avgnucCy5];  

end

    clear sumnucdapi;
    clear avgnucfitc; 
    clear avgnucCy5; 

%The saving of the processed data now happens in the looper
%The whole well is saved as one mat data file
