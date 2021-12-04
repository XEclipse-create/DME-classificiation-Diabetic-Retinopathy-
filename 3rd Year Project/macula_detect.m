function [BW1, BW2, BW3] = macula_detect(I, OD)
I = imread(fullfile('E:\3rd Year Project\images','image013.png'))
OD = imread(fullfile('E:\3rd Year Project\images','image013_OD.png'))
%Resize given image to 576x720
I2=imresize(I, [576 720]); 
%Convert image to HSV
I2_hsv = rgb2hsv(I2);
I2_brightness = I2_hsv(:,:,3);

%Set value of all black pixels (those below given threshold) to 0
I2_black = (I2_brightness<=10);
I2_brightness = double(I2_brightness).*~I2_black;

%Resize optic disc imagee to 576x720
OD2=imresize(OD, [576 720]);

%Set value of all black pixels (those below given threshold) to 0
OD2_black = (OD2<=7);
OD2 = double(OD2).*~OD2_black;

%Check whether optic disc is located in left half or right half of image
od_sum = sum(OD2,1);
od_left = sum(od_sum(1:360));
od_right = sum(od_sum(351:720));

%Determine extreme reach of optic disc, with safety margin
if od_right>od_left
    for i = 361:720
        if od_sum(i)>0
            extreme_left = i-20;
            break;
        end
    end
else
    for i = 1:360
        if od_sum(361-i)>0
            extreme_right = 341-i;
            break;
        end
    end
end

%Determine extreme points of optic disc, without safety margin
for i = 1:720
        if od_sum(i)>0
            od_exleft = i;
            break;
        end
end
for i = 1:720
        if od_sum(721-i)>0
            od_exright = 721-i;
            break;
        end
end

od_sumh = sum(OD2,2);

for i = 1:576
        if od_sumh(i)>0
            od_extop = i;
            break;
        end
end
for i = 1:576
        if od_sumh(577-i)>0
            od_exbottom = 577-i;
            break;
        end
end

%Determine optic disc diameter and centre points
dd = round(((od_exright - od_exleft)+(od_exbottom - od_extop))/2);
od_ycentre = round((od_extop+od_exbottom)/2);
od_xcentre = round((od_exleft+od_exright)/2);
od_centre = [od_xcentre, od_ycentre];

%determine extreme points of fundus image (exclusing black background)
i2_sum = sum(I2_brightness, 1);

if od_right>od_left
    for i = 1:720
        if i2_sum(i)>0
            beg_left = i;
            break;
        end
    end
else
    for i = 1:720
        if i2_sum(721-i)>0
            beg_right = 721-i;
            break;
        end
    end
end

%Display Green Channel of Fundus Image
GreenC=I2(:,:,2);

%Display negative of Green Fundus Image
Ginv2=imcomplement(GreenC); 

%Perform Adaptive Histogram Equalisation on Green Fundus Image
Gadpt_his3=adapthisteq(Ginv2); 

%Perform morphological opening to enhance optical disc
se = strel('ball',8,8);
Gopen4=imopen(Gadpt_his3,se); 

%Subtract Enhanced fundus image to leave behind only the blood vessels
G_Odisk_R5=Gadpt_his3-Gopen4; 


%Perform thresholding
G_BW6 = im2bw(G_Odisk_R5,0.105);

%Remove all obj smaller than pixels value by morphological opening
G_BWareaopen7 = bwareaopen(G_BW6,300);

%Locate Horizontal Edges of Blood Vessels
G_hedge = edge(G_BWareaopen7, 'Sobel', 'horizontal');

%Perform Adaptive Histogram Equalisation on Red Channel of Fundus Image
I2_red = adapthisteq(I2(:,:,1));

%Set value of all black pixels to 0
I2_redzero = (I2_red<=20);

%Locate ROI for macula as defined by equation in paper
if od_left>od_right
    roi_macula = I2_red(od_ycentre-0.5*dd:od_ycentre+2*dd, od_xcentre+dd:od_xcentre+3*dd); 
else
    roi_macula = I2_red(od_ycentre-0.5*dd:od_ycentre+2*dd, od_xcentre-3*dd:od_xcentre-dd);
end

%Define minimum mean and location of centre of macula 
min_mean = inf;
x_centre = 0;
y_centre = 0;


%In given ROI, locate centre of macula by scanning using rectangular window
%of fixed size. Locate the recatngular window which has minimum mean
%brightness. Declare the centre of that rectangle as the centre of the
%macula.
if od_left<od_right
for i=round(od_ycentre-0.5*dd):od_ycentre+2*dd-5
    for j=od_xcentre-3*dd:od_xcentre-dd
        mask = ~I2_redzero(i:i+10, j:j+10);
        rect = I2_red(i:i+10, j:j+10);
        mean_brightness = mean(rect(mask));
        if mean_brightness<min_mean
           min_mean = mean_brightness;
           y_centre = i+5;
           x_centre = j+15;
        end
    end
end
else
    for i=round(od_ycentre-0.5*dd):od_ycentre+2*dd-5
        for j=od_xcentre+dd:od_xcentre+3*dd-5
            mask = ~I2_redzero(i:i+10, j:j+10);
            rect = I2_red(i:i+10, j:j+10);
            mean_brightness = mean(rect(mask));
            if mean_brightness<min_mean
                min_mean = mean_brightness;
                y_centre = i+5;
                x_centre = j+15;
            end
        end
    end
end

%Plot the centre of the macula on top of the Fundus Image
centre = [x_centre, y_centre];

%Display negative of Red Channel Fundus Image 
red_comp = imcomplement(I2_red);

%Perform Adaptive HIstogram Equalisation on Negative Red Channel
red_enhance = adapthisteq(red_comp);

%Calculate Threshold for given image using Otsu's method
level = graythresh(red_enhance);

%Binarise given image using fixed threshold of 0.8
red_thresh = imbinarize(red_enhance,0.8);

%Remove all obj smaller than pixels value by morphological opening
red_areaopen = bwareaopen(red_thresh,300);
%Subtract the previously obtained blood vessels from binary image
red_wobv = red_areaopen - G_BWareaopen7;

%Set all pixel values less than 1 to 0 and display binary image of macula
red_wobv(red_wobv<1)=0;

%Draw obtained prelimiary centre of macula on top of binary image
viscircles(centre, 1, 'EdgeColor', 'r');

%Define new ROI centered on macula using preliminary centre
red_roi = red_wobv(y_centre-dd:y_centre+dd, x_centre-dd:x_centre+dd);

%Locate boundaries of macula in ROI
macula_edges = bwboundaries(red_roi, 'noholes');


%Locate largest continuous boundary and declare it as perimeter of macula
max_length = 0;
for k = 1:length(macula_edges)
   boundary = macula_edges{k};
   if length(boundary) > max_length
       max_length  = length(boundary);
       boundary_large = boundary;
   end
end

%Obtain new Binary Image using boundary of macula as polyshape
[m, n] = size(red_roi);
macula_mask = poly2mask(boundary_large(:,2), boundary_large(:,1), m, n);

%Perform dilation on Binary image to smoothen images
se = strel('disk',5,8);
macula_smooth = imdilate(macula_mask, se);


%Find boundaries to smoothened macula
boundary_smooth = bwboundaries(macula_smooth, 'noholes');

%Compute and display centroid of the smoothened macula
c = regionprops(macula_smooth, 'centroid');
centroid = c.Centroid;

%Compute pixel values of centroid and boundary of macula for original
%fundus image
real_edges = boundary_smooth{1} + [y_centre - dd, x_centre - dd];
im_centroid = centroid + [x_centre - dd, y_centre - dd];

BW1 = createCirclesMask([576,720], im_centroid, 0.33*dd);
BW2 = createCirclesMask([576,720], im_centroid, dd);
BW3 = createCirclesMask([576,720], im_centroid, 2*dd);



