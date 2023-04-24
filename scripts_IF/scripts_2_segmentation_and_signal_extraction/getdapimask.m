function da_ma=getdapimask(da_bs,nucr)

%% threshold mask
da_mas=ThreshImage(da_bs);

%% edge mask
%%
da_bb=[da_bs(:,1+nucr:-1:1+1),da_bs,da_bs(:,end-1:-1:end-nucr)];
da_bb=[da_bb(1+nucr:-1:1+1,:);da_bb;da_bb(end-1:-1:end-nucr,:)];

%%
da_bdx=2*da_bb(1+nucr:end-nucr,1+nucr:end-nucr)-da_bb(1+nucr:end-nucr,1:end-2*nucr)-da_bb(1+nucr:end-nucr,1+2*nucr:end);
da_bdy=2*da_bb(1+nucr:end-nucr,1+nucr:end-nucr)-da_bb(1:end-2*nucr,1+nucr:end-nucr)-da_bb(1+2*nucr:end,1+nucr:end-nucr);
da_bdru=2*da_bb(1+nucr:end-nucr,1+nucr:end-nucr)-da_bb(1:end-2*nucr,1+2*nucr:end)-da_bb(1+2*nucr:end,1:end-2*nucr);
da_bdlu=2*da_bb(1+nucr:end-nucr,1+nucr:end-nucr)-da_bb(1:end-2*nucr,1:end-2*nucr)-da_bb(1+2*nucr:end,1+2*nucr:end);
da_bd=da_bdy+da_bdx+da_bdru+da_bdlu;
% da_bdg=mat2gray(da_bd);
% da_mag0=da_bd>min(da_bd(da_bdg>graythresh(da_bdg)));
da_mag0=ThreshImage(da_bd);
th=0;%th=th0/4;
da_bdl=((da_bdy>th)&(da_bdx>th)&(da_bdru>th)&(da_bdlu>th));
da_mag=da_mag0&da_bdl;

%% combine
da_ma0=imfill((da_mas&da_mag),'holes');

%% segmentation by local maxima
da_bsd=[da_bs(:,1+1),da_bs,da_bs(:,end-1)];
da_bsd=[da_bsd(1+1,:);da_bsd;da_bsd(end-1,:)];

da_bsf=zeros(size(da_bs));
for i=0:2
    for j=0:2
        tempdiff=(da_bs-da_bsd(1+i:end-2+i,1+j:end-2+j));
        da_bsf=da_bsf+(tempdiff>=0);
    end
end
da_pk=(da_bsf>=8)&(da_ma0);
da_pkd=imdilate(da_pk,strel('disk',floor(nucr*2/3),0));
da_wa=watershed(bwdist(da_pkd))>0;
da_ma1=(da_ma0&da_wa);

%% segmentation of da_mas
% da_wa2=watershed(bwdist(da_ma1))>0;
% da_ma2=(da_mas&da_wa2);

%% filtering
x=0.5;  %change this to change the small debris cutoff.  default is x=1 
da_ma=bwareaopen(da_ma1,round(x*(nucr.^2)));  %if nucr=8, code will remove any cells w are <64 pixels
% figure
% imshow(da_ma)
% title('NUCLEAR MASK')
