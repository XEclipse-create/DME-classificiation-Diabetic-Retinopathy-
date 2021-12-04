
clc; clear; close all;
I = imread(fullfile('E:\3rd Year Project\images','21_training.tif'));

%% Preprossecing

red = I(:,:,1);
gaussianFilter = fspecial('gaussian',35, 25);
gaus = imfilter(red, gaussianFilter,'symmetric');

levels = multithresh(red, 10);
BW = imquantize(gaus,levels);
BW(BW<11)=0;

se = strel('disk',2);
ED = imopen(BW, se);

se = strel('disk',13);
ED2 = imdilate(ED, se);

img_edges = edge(ED2, 'Canny');

se = strel('disk',1);
ED3 = imdilate(img_edges,se);

%% Circle Hough Transform

radius_range = [40, 50];
r_min = radius_range(1);
r_max = radius_range(2);
r_num = 20;
numpeaks = 6;
centers = zeros(r_num * numpeaks, 2);
radii = zeros(size(centers,1),1);
row_num = 0;
for radius = linspace(r_min, r_max, 5)
    % Compute Hough accumulator array for finding circles.
    H = zeros(size(ED3));
    for x = 1 : size(ED3, 2)
        for y = 1 : size(ED3, 1)
            if (ED3(y,x))
                for theta = linspace(0, 2 * pi, 360)
                    a = round(x + radius * cos(theta));                
                    b = round(y + radius * sin(theta));
                    if (a > 0 && a <= size(H, 2) && b > 0 && b <= size(H,1))
                        H(b,a) = H(b,a) + 1;
                    end
                end
            end
        end
    end
    
    % Find peaks in a Hough accumulator array
    threshold = 0.8 * max(H(:)); % which values of H are considered to be peaks
    nHoodSize = floor(size(H) / 100.0) * 2 + 1; % Size of the suppression neighborhood, [M N]
    peaks = zeros(numpeaks, 2);
    num = 0;
    while(num < numpeaks)
        maxH = max(H(:));
        if (maxH >= threshold)
            num = num + 1;
            [ra,c] = find(H == maxH);
            peaks(num,:) = [ra(1),c(1)];
            rStart = max(1, ra - (nHoodSize(1) - 1) / 2);
            rEnd = min(size(H,1), ra + (nHoodSize(1) - 1) / 2);
            cStart = max(1, c - (nHoodSize(2) - 1) / 2);
            cEnd = min(size(H,2), c + (nHoodSize(2) - 1) / 2);
            for i = rStart : rEnd
                for j = cStart : cEnd
                        H(i,j) = 0;
                end
            end
        else
            break;          
        end
    end
    peaks = peaks(1:num, :);        

    if (size(peaks,1) > 0)
        row_num_new = row_num + size(peaks,1);
        centers(row_num + 1:row_num_new,:) = peaks;
        radii(row_num + 1:row_num_new) = radius;
        row_num = row_num_new;       
    end
end

centers = centers(1:row_num,:);
radii = radii(1:row_num);


summ = sum(ED2(:));
if summ > 80000
    i = 5;
else
    i = 10;
end

%% Draw this circle on the image

figure();
imshow(I);
hold on;
rmax = 1;
r = radii(i);
center_x = centers(i, 2);
center_y = centers(i, 1);
theta = linspace(0, 2 * pi, 360);
xx = center_x + r * cos(theta);
yy = center_y + r * sin(theta);
plot(xx, yy,'r', 'LineWidth', 1);
hold off;
[r center_y center_x]