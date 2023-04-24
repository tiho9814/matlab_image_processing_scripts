function exportND2toTIF_Jian(nd2File,dirToSave)

% EXPORTND2TOTIF  Export ND2 file to TIF files which work with Mingyu's code
% Naming convention: row_col_1_Pos5 ChannelName_Frame
% edited by Jian 170525

if ~exist(dirToSave,'dir')
    mkdir(dirToSave)
end

%% Resolve the well column and row from the nd2 filename
[~,fn] = fileparts(nd2File);

% Find a capital letter, followed by two digits
startIdx = regexp(fn,'[A-Z][0-9][0-9]');

% Convert the first letter into a number
wellRow = int8(fn(startIdx)) - 64; %'A' = 65

% Convert the next two digits into the column number
wellCol = str2double(fn(startIdx+1:startIdx+2));

%% Export the images
% nd2r = bfGetReader(nd2File);
bfReaderJian = BioformatsImage(nd2File);


%% number frames and channels
% numFrames = nd2r.getSizeT;
% numChans  = nd2r.getSizeC;
% numSeries = nd2r.getSeriesCount;

% Get the channel name for the metadata and define the new names for the
% exported images
% md = nd2r.getMetadataStore;
% channelFNPrefix = cell(1,numChans);
% for iC = 1:numChans
%     channelFNPrefix{iC} = char(md.getChannelName(0,iC - 1));
% end

%Convert Well Name to number
% hWB = waitbar(0,'Converting...');
% AH for iFrame = 1:numFrames
% AH waitbar(iFrame/numFrames, hWB);

for xyLoc = 1:bfReaderJian.seriesCount     % SITE

   for iC = 1:bfReaderJian.sizeC           % Channel
        
        for iFrame = 1:bfReaderJian.sizeT  % Frame
            
            img = bfReaderJian.getXYplane(iC,xyLoc,iFrame);
            
        % set the current file name
            currFN = sprintf('%d_%d_%d_%s_%d.tif',wellRow, wellCol, xyLoc, bfReaderJian.channelNames{iC},iFrame);

        % write the data to the tiff file
            imwrite(img,fullfile(dirToSave,currFN),'tif','Compression','None');
        end
        
    end
    
end


end