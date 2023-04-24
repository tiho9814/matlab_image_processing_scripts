function rpoint = rpoint_check_V2(dhb)

% identify R point
% Mingwei Min 5.15.2017
% 
% version 2: use peak of CDK2 slope as primary identification
% test cdk2 level(+- a timewindow) of cadidates 
% refine around it
% Mingwei Min 5.17.2017

rpoint = 0;
if sum(~isnan(dhb))>20
    xrange          = find(~isnan(dhb));
    dhb_raw         = dhb(xrange);
    dhb_smoothed    = smooth(dhb_raw,20,'loess');

    abnormal_test   = max(dhb_smoothed)<3; %& max(diff(dhb_raw))<0.2; % test if the CDK2 activity exceeds normal range
    
    if abnormal_test
        gate            = zeros(1,length(xrange));
        slopewindow     = 10;
        dhb_slope       = getslope_forward_avg(dhb_smoothed,4:8);
        rises           = peakdet(dhb_slope,0.0001);
        
        % go through each reflection point
        for j = 1:size(rises,1)
            
            % identify R point area
            if rises(j,1)>10 & rises(j,1)<=length(dhb_smoothed)-2*slopewindow
                testpoint       = rises(j,1);
                pastheight      = mean(dhb_smoothed(max(testpoint-slopewindow,1):testpoint))<0.7;%was 1.5
                futureheight    = max(dhb_smoothed(testpoint:testpoint+slopewindow))>dhb_smoothed(testpoint)+0.15; %0.15
                futureheight2   = mean(dhb_smoothed((testpoint+slopewindow):(testpoint+2*slopewindow)))>0.6;
                nextlow         = find(dhb_smoothed((testpoint+slopewindow):end)<0.6,1,'first');
                if nextlow
                    futurelow   = nextlow>30;
                    betweenlow  = max(dhb_smoothed(testpoint:testpoint+nextlow))>1.1;
                    not_spike   = futurelow & betweenlow;
                else
                    not_spike   = 1;
                end
                
                % refine R point
                while pastheight & futureheight & futureheight2 & double(not_spike) & testpoint>=10
                    rpoint          = testpoint;
                    testpoint       = testpoint-1;
                    pastheight      = mean(dhb_smoothed(max(testpoint-slopewindow,1):testpoint))<0.7;%was 1.5
                    futureheight    = max(dhb_smoothed(testpoint:testpoint+slopewindow))>dhb_smoothed(testpoint)+0.15; %0.15
                    futureheight2   = mean(dhb_smoothed((testpoint+slopewindow):(testpoint+2*slopewindow)))>0.6;
                    nextlow         = find(dhb_smoothed((testpoint+slopewindow/2):end)<0.6,1,'first');
                    if nextlow
                        futurelow   = nextlow>30;
                        betweenlow  = max(dhb_smoothed(testpoint:testpoint+nextlow))>1.1;
                        not_spike   = futurelow & betweenlow;
                    else
                        not_spike   = 1;
                    end
                end
                
                if rpoint
                    rpoint = xrange(rpoint);
                    break;
                end
            end
        end   
    end
end


