function imshowc(varargin)

I(:,:,1)=im2double(varargin{1})*1;
I(:,:,2)=zeros(size(varargin{1}));
I(:,:,3)=zeros(size(varargin{1}));

if (length(varargin)>2)
    I(:,:,2)=im2double(varargin{2})*1;
end

if (length(varargin)>3)
    I(:,:,3)=im2double(varargin{3})*1;
end

if (varargin{end}==1)
    for i=1:3
        I(:,:,i)=imadjust(mat2gray(I(:,:,i)));
    end
end
imshow(I);