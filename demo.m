addpath('./core/');

filename = 'inputs/1.jpg'; 
% parameters
Tac = 165;
Tni = 0.5;

% read image
disp('read image------------------------------------------------');
I = imread(filename);
figure;imshow(I);
% circle detection
disp('circle detetion-------------------------------------------');
[circles, ~,~] = circleDetectionByArcsupportLS(I, Tac, Tni);
% display
disp('show------------------------------------------------------');
circles
disp(['number of circles��',num2str(size(circles,1))]);
disp('draw circles----------------------------------------------');
dispImg = drawCircle(I,circles(:,1:2),circles(:,3));
figure;
imshow(dispImg);
