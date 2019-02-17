% corner detection
% eigen value points to the directions of changes; bigger eigen value, bigger change

% load image
Img = imread('sfu.jpg');
I = rgb2gray(Img);
I = padarray(I,[1 1], 0, 'both');

% compute gradient
% gx = [1, 0, -1; 2, 0, -1; 1, 0, -1];
gx = [-1, 0, 1; -2, 0, 2; -1, 0, 1];
gy = gx';
Ix = conv2(I, gx, 'same');
Iy = conv2(I, gy, 'same');

% apply gaussian smooth
Ix2 = Ix .^ 2;
Iy2 = Iy .^ 2;
IxIy = Ix .* Iy;

Ix2 = imgaussfilt(Ix2);
Iy2 = imgaussfilt(Iy2);
IxIy = imgaussfilt(IxIy);


[row, col] = size(I);
L1 = zeros(row, col);
L2 = zeros(row, col);
Ix2p = padarray(Ix2, [1, 1]);
Iy2p = padarray(Iy2, [1, 1]);
Ixyp = padarray(IxIy, [1, 1]);

M11 = Ix2p;
M12 = Ixyp;
M21 = Ixyp;
M22 = Iy2p;

sumFilter = [1 1 1; 1 1 1; 1 1 1];

M11sum = conv2(M11, sumFilter, 'same');
M12sum = conv2(M12, sumFilter, 'same');
M21sum = conv2(M21, sumFilter, 'same');
M22sum = conv2(M22, sumFilter, 'same');

E = (1/2) * (M11sum + M22sum) - sqrt(4 * (M12sum + M21sum) + (M11sum - M22sum) .^ 2);
% thres = 0.5;

% [row, col] = size(I);
% thres = 50000;
% 
% E = zeros(row, col); 
% for i = 2 : 1 : row - 1 
% 	for j = 2 : 1 : col - 1
%         Ix2_sum = sum(sum(Ix2(i-1 : i+1, j-1 : j+1)));
%         Iy2_sum = sum(sum(Iy2(i-1 : i+1, j-1 : j+1)));
%         IxIy_sum = sum(sum(IxIy(i-1 : i+1, j-1 : j+1)));
%         % build m for pixel
%         M = [Ix2_sum, IxIy_sum; IxIy_sum, Iy2_sum];
%         % calculate eigenvalues
% %         E(i,j) = min(eig(M));
%         % threshold lambda
%         lambda_pos = (1/2) * (M(1, 1) + M(2 + 2)) + sqrt(4 * (M(1, 2) + M(2, 1)) + (M(1, 1) - M(2, 2)) ^ 2);
%         lambda_neg = (1/2) * (M(1, 1) + M(2 + 2)) - sqrt(4 * (M(1, 2) + M(2, 1)) + (M(1, 1) - M(2, 2)) ^ 2);
%         if (lambda_pos > thres && lambda_neg > thres)
%             E(i, j) = min(lambda_pos, lambda_neg);
%         end
% 	end
% end

% disp(E);
imshow(E > 50000);

