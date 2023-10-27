function [ outputImg ] = drawCircle( Img, Centers, Rs )
%draw circles in the Img
%Img        : Img is 1 or 3 channels image
% Centers is: N x 2,(col,row),namely(x,y)
% Rs      is: N    

%{
%method 3
%show circles
h = size(Img,1);
w = size(Img,2);
if(size(Img,3) == 1) %如果输入图像为单通道图像
outputImg = zeros(h,w,3);
outputImg(:,:,1) = Img;
outputImg(:,:,2) = Img;
outputImg(:,:,3) = Img;
else if (size(Img,3) == 3)%如果输入图像为3通道图像
        outputImg = Img;
    end
end
N     = size(Centers,1);
for i = 1 : N
    centerx = floor(Centers(i,2));
    centery = floor(Centers(i,1));
    for j = 1:360        
        x = centerx+floor(Rs(i).*cos(pi./180.*j));%以图像高为x轴       
        y = centery+floor(Rs(i).*sin(pi./180.*j));%以图像宽为y轴
        if( x>=1 && x<=h && y>=1 && y<= w )
            outputImg(x,y,1) = 255;
            outputImg(x,y,2) = 0;
            outputImg(x,y,3) = 0;
        end
    end
    if( centerx>=1 && centerx<=h && centery>=1 && centery<= w )
     outputImg(centerx,centery,1) = 255;
     outputImg(centerx,centery,2) = 0;
     outputImg(centerx,centery,3) = 0;
    end
end
outputImg = uint8(outputImg);
%}

%method 1
%show circles
if (size(Img,3) == 1)
outputImg = zeros(size(Img,1),size(Img,2),3);
outputImg(:,:,1)=Img;
outputImg(:,:,2)=Img;
outputImg(:,:,3)=Img;
else
    outputImg = Img;
end
[y, x] = ndgrid(1 : size(Img, 1), 1 : size(Img, 2));
x = x(:); y = y(:);
for i = 1 : size(Centers, 1)
    idx = abs(sqrt((x - Centers(i, 1)) .^ 2 + (y - Centers(i, 2)) .^ 2) - Rs(i)) <= 1;
    outputImg(sub2ind(size(outputImg), y(idx), x(idx), ones(sum(idx), 1))) = 0;
    outputImg(sub2ind(size(outputImg), y(idx), x(idx), 2 * ones(sum(idx), 1))) = 255;
    outputImg(sub2ind(size(outputImg), y(idx), x(idx), 3 * ones(sum(idx), 1))) = 0;
    %draw circular centers golden 235 199 16
    outputImg(round(Centers(i,2)+0.5),round(Centers(i,1)+0.5),:)=[255 255 255];
end
outputImg = uint8(outputImg);


%{
%method 2
%show circles
figure(); 
hold on;
imshow(Img);
for i = 1 : size(Centers, 1)
    viscircles(Centers,  Rs, 'EdgeColor', 'g', 'DrawBackgroundCircle', false, 'LineWidth', 1);
end
%}
end