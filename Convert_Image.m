function [Image2] = Convert_Image(Image,H,W)

% H=H+1;
% W=W+1;

Image2=zeros(H,W);

Imax=max(max(Image));

Imin=min(min(Image));


Image2=round(((Image-Imin)/(Imax-Imin))*255);
% Image2=flipud(Image2);
% Image2=fliplr(Image2);

% for i=1:H
%     for j=1:W
%         
%         I=Image(i,j);
%         Image2(i,j)=((I-Imin)/(Imax-Imin))*255;
%         Image2(i,j)=round(Image2(i,j));
%     end
% end
end




