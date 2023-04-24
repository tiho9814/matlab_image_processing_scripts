function cells = cellinfo_V3_DHB_mCherry(datadir, row, col, site, nuc_col, cyto_col, framerate)

% version 2, added in 
% 1. test for mitosis skip (potential senescence entry)
% 2. and identification of R point
% 3. cdk2 state of mother
% 4. test if the trace is noisy
% Mingwei Min, 2017.5.17

% addition in 2017.5.23: set NA for mother quality instead of empty

% addition in 2017.8.23: handle empty trace

% Version 3,
% improved cdk2 state categorization. 2018.5.21
% added in filter to remove daugther cell that lose the sensor

% addition in 2018.8.13: improve detection of emerging cells - if a cell
% has a noisy cdk2 increase followed by a real increase, it would be
% classified as a low cell before. Now it's an inc cell.

shot        = [num2str(row),'_', num2str(col), '_', num2str(site)];
load([datadir,'tracedata_',shot,'.mat'],'tracedata','genealogy','jitters');

cellnum     = length(genealogy);
cells       = struct([]);
callthresh  = 0.5;
callthresh2 = 0.7;

for i=1:cellnum
    cells(i).raw    = tracedata(i,:,:);
    cells(i).raw    = squeeze(cells(i).raw);
    raw             = cells(i).raw;
    trace_length    = sum(~isnan(raw(:,1)));
    if trace_length>0
        trace_start     = find(~isnan(raw(:,1)),1,'first');
        trace_end       = find(~isnan(raw(:,1)),1,'last');
    else
        trace_start     = 0;
        trace_end       = 0;
    end
    
    cells(i).cdk2_trace     = raw(:,cyto_col)./raw(:,nuc_col);
    cdk2_trace      = cells(i).cdk2_trace;
    if trace_length>0
        cdk2_trace(trace_start:trace_end)   = smooth(cdk2_trace(trace_start:trace_end));
    end
    
    %% set default sensor quality
    cells(i).momQuality = 2; % 2 as the sensor quality cannot be decided    
    
    %% test if the CDK2 sensor is lost
    if (nanmean(raw(:,cyto_col))+nanmean(raw(:,nuc_col)))>20
        cells(i).sensor  = 1;
    else
        cells(i).sensor  = 0;
    end
    
    %% test if the cdk2 trace is noisy or beyond normal range
    if max(cdk2_trace)>3 | min(cdk2_trace)<0.0001 | max(diff(cdk2_trace))>0.2
        cells(i).noisy  = 1;
    else
        cells(i).noisy  = 0;
    end
    
    %% if the cell is a mother
    if ismember(i,genealogy)
        cells(i).last_mito  = trace_end+1;
        cells(i).daughter   = find(genealogy==i);
        
        % if it's a mother, test if the sensor is working fine
        mom_quality_check   = min(trace_start, max(1,trace_end - framerate*3));
        if sum(~isnan(cdk2_trace(mom_quality_check:end)))>5
            cells(i).momQuality = max(cdk2_trace(mom_quality_check:end))>1;
        end
        
        
        % if it's a mother, it does not skip mitosis
        cells(i).mito_skip  = 0;
    else
        cells(i).last_mito  = 0;
        cells(i).daughter   = 0;
        % if it's not a monther, check if it goes from
        % G2-like to G1-like CDK2 state without a mitosis
        if cells(i).sensor
            cells(i).mito_skip = mitosis_skip_check(cdk2_trace, raw(:,nuc_col), raw(:,3),framerate,~isnan(genealogy(i)));
        else
            cells(i).mito_skip = 0;
        end
        
    end

    %% if the cell is a daughter
    if ~isnan(genealogy(i))
        cells(i).first_mito  = trace_start;
        cells(i).mother      = genealogy(i);
        cells(i).mom_cdk2    = cells(genealogy(i)).cdk_state;
    % catergorize CDK2 group  
        if cells(genealogy(i)).momQuality==0 % if mother sensor fail
           categorization   = 3; % do not catergorize
           cells(i).sensor  = 0; % its own sensor fail too
        else
            calltime    = trace_start + 3*framerate;
            if calltime>trace_end
                categorization  = 3;
            elseif max(cdk2_trace(trace_start+framerate:trace_start+2*framerate))>0.8 %if cdk2 activity is above 0.8 right after mitosis
                categorization      = 3; % do not catergorize
                cells(i).sensor     = 0; % its sensor fail
            else            
                categorization  = (cdk2_trace(calltime))>callthresh;
                if categorization == 1 
                    test_point_min = min(calltime+4*framerate,trace_end);
                    if cdk2_trace(test_point_min)<callthresh2
                        categorization = 3;
                    end
                else
                    if max(cdk2_trace(calltime:trace_end))>0.6
                        frames_to_check = calltime:trace_end; % test if the increase of CDK2 activity is noise
                        inc_frame = find(cdk2_trace(frames_to_check)>0.6);
                        test_frame = inc_frame(round(length(inc_frame)/2));
                        if length(inc_frame)/2*(trace_end - frames_to_check(test_frame))> 0.8 %if most of frames after the first >0.6 frame have CDK2 activity>0.6
                            categorization = 2;
                        end
                    end
                end
            end
        end
    else
        cells(i).first_mito = 0;
        cells(i).mother     = 0;
        categorization      = 3;
    end

    %% test for R point
    if categorization>0
        cells(i).rpoint = rpoint_check_V2(cdk2_trace);
    else
        cells(i).rpoint = 0;
    end
    
    cells(i).trace_length   = trace_length;
    cells(i).trace_start    = trace_start;
    cells(i).trace_end      = trace_end;
    cells(i).cdk_state      = categorization;
end
save([datadir,'cellinfo_',shot,'.mat'],'cells');
