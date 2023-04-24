function [ImageMask,th,bg]=ThreshImage(ImageOld)
%ngps:number of groups
%cutafter:cut after which group (from left)

%% make judgement
TempSeries=ImageOld(:);

%% get the histogram of the image values
tmax=max(TempSeries);
tmin=min(TempSeries);
nbin=200;tbin=(tmax-tmin)/nbin;
tmin=tmin+tbin/2;tmax=tmax-tbin/2;
[n,xout]=ksdensity(TempSeries,tmin:tbin:tmax);
% added to look at the histogram
% figure
% plot(xout,n)
% remove later
gp=max([2,ceil(nbin/50)]);
ng=getcurvature(n,gp);

%% get background value
negpeak=regionprops(bwlabel(ng<0),'PixelIdxList','Area');
negth=prctile(ng,25)-1.5*iqr(ng);
Ibg00=zeros(1,size(negpeak,1));AC=Ibg00;DD=Ibg00;
for cc=1:size(negpeak,1)
    AC(cc)=negpeak(cc).Area;
    tempeakvalues=ng(negpeak(cc).PixelIdxList);
    [DD(cc),I]=min(tempeakvalues);
    Ibg00(cc)=negpeak(cc).PixelIdxList(I);
end
neglo=((DD<=negth)|(AC>=gp));Ibg00=Ibg00(neglo);DD=DD(neglo);
if isempty(Ibg00)
    [dummy,Ibg00]=min(ng);
end
Ibg00=[Ibg00,size(ng,2)-gp];
[dummy,I_DD]=min(DD);
Ibg0=Ibg00(I_DD);
bg0=xout(Ibg0);%approximate value
TempSeries2=TempSeries((TempSeries>(bg0-5*gp*tbin))&(TempSeries<(bg0+5*gp*tbin)));
[n_bg,xout_bg]=ksdensity(TempSeries2);
[dummy,Ibg]=max(n_bg);%better value
xfit=xout_bg(Ibg-1:Ibg+1)';yfit=n_bg(Ibg-1:Ibg+1)';
det_n=det([xfit.^2,xfit,[1;1;1]]);
a=det([yfit,xfit,[1;1;1]])/det_n;
b=det([xfit.^2,yfit,[1;1;1]])/det_n;
bg=-b/a/2;%best value

%% get threshold
TempSeries3=TempSeries((TempSeries>(bg-gp*tbin))&(TempSeries<(xout(Ibg00(I_DD+1))+10*gp*tbin)));
[n_th,xout_th]=ksdensity(TempSeries3,bg:(xout(Ibg00(I_DD+1))+10*gp*tbin-bg)/100:(xout(Ibg00(I_DD+1))+10*gp*tbin));
[dummy,Ibg]=min(abs(xout_th-bg));
ng_th=getcurvature(n_th,gp);
partng=ng_th(Ibg:end);
pospeak=regionprops(bwlabel(partng>0),'PixelIdxList','Centroid','Area');
Ith00=zeros(1,size(pospeak,1));AC=Ith00;
for cc=1:size(pospeak,1)
    AC(cc)=sum(partng(pospeak(cc).PixelIdxList));
    %AC(cc)=pospeak(cc).Area;
    Ith00(cc)=round(pospeak(cc).Centroid(1));
    %Ith00(cc)=(sum((pospeak(cc).PixelIdxList)'.*(partng(pospeak(cc).PixelIdxList)))/AC(cc));
end
[dummy,I]=max(AC);Ith0=floor(Ith00(I));Ithr=(Ith00(I)-Ith0);
s1=xout_th(Ith0+Ibg-1);s2=xout_th(Ith0+Ibg);
th=s1+Ithr*(s2-s1);

%% get mask
ImageMask=single(ImageOld>th);