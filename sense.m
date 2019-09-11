function im_full = sense( ima, map, RX, RY )
%
% im = sense( ima, map, RX, RY ) with caipirinha
%
%  ima -- aliased images, undersampled at Rx and Ry, at full FOV
%  map -- maps of the coil sensitivities
%  Rx  -- x acceleration factor
%  Ry  -- y acceleration factor
%
%  im  -- SENSE reconstruction of ima
%
%  this assumes that the noise correlation matrix is I
%

%
% written by John Pauly
% (c) Board of Trustees, Leland Stanford Jr University, 2011
%
ima_dum=zeros(size(map));

[NX NY L] = size(map );

NRX = NX/RX;
NRY = NY/RY;
%%%%%%%%%%%%%%%%%%%%%%%%
% estimate sensitivites
%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Estimating Sensitivities...\n');
sosimg = sqrt(sum(abs(map.^2),3));
map= (map./repmat(sosimg,[1 1 L]));
im = zeros(NX,NY);
ima_full = zeros(NX,NY);
shift=0;
% for ll=1:L
%     ima_temp=fft2c(ima(:,:,ll));
% %     map_temp=fft2c(map(:,:,ll));
% %     map_temp=circshift(map_temp,[0,shift]);
% %     ima_temp=circshift(ima_temp,[0,shift]);
%     ima_dum(1:RX:NX,1:RY:NY,ll)=(ima_temp);
%     ima_full(:,:,ll)=ifft2c(ima_dum(:,:,ll));
% %     map(:,:,ll)=ifft2c(map_temp);
% end
ima_full=ima;
for ii=1:NX,
    for jj=1:NY,
        if abs(map(ii,jj,1)) <= 0,
            im(ii,jj) = 0;
        else
            for LX=0:RX-1,
                for LY=0:RY-1,
                    ndx = mod((ii-1)+LX*NRX,NX)+1;
                    ndy = mod((jj-1)+LY*NRY,NY)+1;
                    CT = map(ndx, ndy, :);
                    CT = CT(:);
                    if ((LX==0) & (LY==0)),
                            s = CT;
                    else if abs(CT(1)) > 1e-6,
                            s = [s CT];
                        end;
                    end;
                end;
            end;
            scs = s'*s;
            scsi = inv(scs);
            m = ima_full(ii,jj,:);
            m = m(:);
            mr = scsi*s'*m;
%             mr = pinv(s)*m;
            im(ii,jj) = mr(1);
            
        end;
    end;
end
im_shift=fft2c(im);
im_shift=circshift(im_shift,[0,-shift]);
im_full=ifft2c(im_shift);
