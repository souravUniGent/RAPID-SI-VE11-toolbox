function I=image_smooth(I,sigma)

siz=2*sigma;%sigma;
% Make 1D gaussian kernel
range=-ceil(siz/2):ceil(siz/2);
H = exp(-(range.^2/(2*sigma^2)));
H = H/sum(H(:));
% Filter each dimension with the 1D gaussian kernels\
if(ndims(I)==1)
    I=imfilter(I,H', 'same' ,'replicate');
elseif(ndims(I)==2)
    Hx=reshape(H,[length(H) 1]);
    Hy=reshape(H,[1 length(H)]);
    I=imfilter(imfilter(I,Hx, 'same' ,'replicate'),Hy, 'same' ,'replicate');
elseif(ndims(I)==3)
    if(size(I,3)<4) % Detect if 3D or color image
        Hx=reshape(H,[length(H) 1]);
        Hy=reshape(H,[1 length(H)]);
        for k=1:size(I,3)
            I(:,:,k)=imfilter(imfilter(I(:,:,k),Hx, 'same' ,'replicate'),Hy, 'same' ,'replicate');
        end
    else
        Hx=reshape(H,[length(H) 1 1]);
        Hy=reshape(H,[1 length(H) 1]);
        Hz=reshape(H,[1 1 length(H)]);
        I=imfilter(imfilter(imfilter(I,Hx, 'same' ,'replicate'),Hy, 'same' ,'replicate'),Hz, 'same' ,'replicate');
    end
end