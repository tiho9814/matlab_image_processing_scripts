function senescent = mitosis_skip_check(dhb, nuc, nuc_h2b, framerate, is_daughter)

% check if the cell goes from a CDK2 high state to a CDK2 low state without
% a mitosis.
% Mingwei Min 5.13.2017

% for i = 123:length(cells)
%     dhb = cells(i).cdk2_trace;
%     nuc = tracedata(i,:,11);
%     nuc_h2b = tracedata(i,:,3);
%     is_daughter = ~isempty(cells(i).mother);
senescent = 0;
if sum(~isnan(dhb))>10

    xrange          = find(~isnan(dhb));
    dhb_raw         = dhb(xrange);
    dhb_smoothed    = smooth(dhb_raw,10,'loess');

    dhb_slope       = getslope_forward_avg(dhb_smoothed,4:8);
    dhb_slope       = dhb_slope*(-1);
    drops           = peakdet(dhb_slope,0.03);
    abnormal_test   = max(dhb_smoothed)<6; % test if the CDK2 activity exceeds normal range
    
    if ~isempty(drops) & abnormal_test
        max_drop        = sortrows(drops,-2);       
        nuc_raw         = smooth(nuc(xrange),10,'loess');
        nucdhb_slope    = getslope_forward_avg(nuc_raw,4:8);
        
        mean_nuc_area       = nanmean(nuc_h2b);
        normalized_nuc_area = nuc_h2b/mean_nuc_area;    
        nuc_area_smoothed   = smooth(normalized_nuc_area,10,'loess'); % smooth normalized 
        nuc_area_slope      = getslope_forward_avg(nuc_area_smoothed(xrange),4:8);
        nuc_area_slope      = nuc_area_slope*(-1);
        max_areaChange      = peakdet(nuc_area_slope,0.02);
        
        % loop through each CDK2 drop
        for j = 1:size(max_drop,1)
            testpoint           = max_drop(j,1);     
            G2test              = max(dhb_smoothed(max(testpoint-2*framerate,1):testpoint))>1 & mean(dhb_smoothed(max(testpoint-2*framerate,1):testpoint))>0.8;
            lengthtest          = length(dhb_smoothed)>(testpoint+2*framerate);
            G1test              = mean(dhb_smoothed((testpoint+2*framerate):min(testpoint+6*framerate,length(dhb_smoothed))))<0.9;
            translocate_test    = mean(nucdhb_slope(testpoint:testpoint+framerate))>0; % excluding the drops from moving into/away from a bright cell. Test for change of nuc intensity. Change of cyto intensity can be too subtle or too noisy to
            mitosisdrop_test    = ~(is_daughter & testpoint<4*framerate); % if it's a daughter and the point is within 20 frames of birth, exclude it - just drop from mitosis
            
            if G2test & lengthtest & G1test & translocate_test & mitosisdrop_test
                senescent       = xrange(testpoint);
                % exclude those did mitosis but mis-tracked (condensed chromatin)
                for k = 1:size(max_areaChange,1)
                    if abs(max_areaChange(k,1)-testpoint)<6
                        senescent   = 0;
                        break;
                    end
                end
                break;
            end
        end
 
% % for testing only
%                 figure, hold on
%                 subplot(3,1,1)
%                 plot(xrange,dhb_raw);
%                 hold on;
%                 plot(xrange,dhb_smoothed);
%                 subplot(3,1,2)
%                 plot(xrange,dhb_slope);
%                 subplot(3,1,3)
%                 plot(xrange,nuc_mean_slope);
%     end
    end

end
