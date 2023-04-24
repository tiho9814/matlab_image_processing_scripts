function ciplot_pro(input,x,colour);
     
% ciplot(lower,upper)       
% ciplot(lower,upper,x)
% ciplot(lower,upper,x,colour)
%
% Plots a shaded region on a graph between specified lower and upper confidence intervals (L and U).
% l and u must be vectors of the same length.
% Uses the 'fill' function, not 'area'. Therefore multiple shaded plots
% can be overlayed without a problem. Make them transparent for total visibility.
% x data can be specified, otherwise plots against index values.
% colour can be specified (eg 'k'). Defaults to blue.

% Raymond Reynolds 24/11/06
% modified by Mingwei Min 8/22/16
% modified again 4/21/16 to handle the whold matric
% mean_measure = nanmean(input,1);
mean_measure = nanmedian(input,1); %modified by MA to plot median instead of mean
std_measure = nanstd(input,1);
input = input(sum(~isnan(input),2)>0,:);
n = size(input,1);

%x = 1:size(input,2);

gate=ones(length(x),1); %remove all NaNs
gate=(~isnan(mean_measure) & (~isnan(x)));
x=x(gate);
mean_measure=mean_measure(gate);
std_measure=std_measure(gate);

lower=mean_measure-2*std_measure./sqrt(n);
upper=mean_measure+2*std_measure./sqrt(n);
if length(lower)~=length(upper)
    error('lower and upper vectors must be same length')
end

if nargin<3
    colour='b';
end

if nargin<3
    x=1:length(lower);
end

% convert to row vectors so fliplr can work
if find(size(x)==(max(size(x))))<2
x=x'; end
if find(size(lower)==(max(size(lower))))<2
lower=lower'; end
if find(size(upper)==(max(size(upper))))<2
upper=upper'; end

line(x,mean_measure,'color',colour);
hold on
h=fill([x fliplr(x)],[upper fliplr(lower)],colour);
set(h,'facealpha',.5,'EdgeColor','none');
