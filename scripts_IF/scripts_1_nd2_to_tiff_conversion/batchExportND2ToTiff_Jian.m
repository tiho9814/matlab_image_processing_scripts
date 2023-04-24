%install BioformatsImage.mltbx before beginning

clearvars
clc
opengl hardware;

dirToProcess ='/Users/tim/Desktop/scripts_IF/example_process/';   % directory containing nd2 files
dirToExport  = '/Users/tim/Desktop/scripts_IF/example_process/output_1_TIFF_files/';  % directory to export individual TIFF files

%p = parpool(2);
if ~exist(dirToExport,'dir')
    mkdir(dirToExport)
end

fileList = dir(fullfile(dirToProcess,'*.nd2'));
fileList = {fileList.name};
list = [1];                  %[1:96] for a full 96-well plate imaged
fileList = fileList(list);
parfor iFile = 1:numel(fileList)  
    
    currFile = fullfile(dirToProcess,fileList{(iFile)});
    disp(currFile)
    
    exportND2toTIF_Jian(currFile,dirToExport); 
end
%close(p)
