%% code based https://github.com/MIT-SPARK/CertifiablyRobustPerception/
%% quasar_stitching_gen_data
load(dataFileName)
R = R_stride';
H = K*R*Kinv;
H = H';
tforms1(2) = projective2d(eye(3));
tforms1(2) = projective2d(H);
xlim=zeros(2,2);
ylim=zeros(2,2);
for i = 1:numel(tforms1)
    [xlim(i,:), ylim(i,:)] = outputLimits(tforms1(i), [1 imageSize(2)], [1 imageSize(1)]);
end
maxImageSize = imageSize;

% Find the minimum and maximum output limits
xMin = min([1; xlim(:)]);
xMax = max([maxImageSize(2); xlim(:)]);

yMin = min([1; ylim(:)]);
yMax = max([maxImageSize(1); ylim(:)]);

% Width and height of panorama.
width  = round(xMax - xMin);
height = round(yMax - yMin);

% Initialize the "empty" panorama.
panorama = zeros([height width 3], 'like', image1);

blender = vision.AlphaBlender('Operation', 'Binary mask', ...
    'MaskSource', 'Input port');

% Create a 2-D spatial reference object defining the size of the panorama.
xLimits = [xMin xMax];
yLimits = [yMin yMax];
panoramaView = imref2d([height width], xLimits, yLimits);
numImages=2;
imgCell={image1Col,image2Col};
for i = 1:numImages
    
    I = imgCell{i};
    
    % Transform I into the panorama.
    warpedImage = imwarp(I, tforms1(i), 'OutputView', panoramaView);
    
    % Generate a binary mask.
    mask = imwarp(true(size(I,1),size(I,2)), tforms1(i), 'OutputView', panoramaView);
    
    % Overlay the warpedImage onto the panorama.
    panorama = step(blender, panorama, warpedImage, mask);
end

figure;
imshow(panorama)
shg

%% plot inlier and outlier correspondences
figure;
stackedImage = cat(2, image1, image2); % Places the two images side by side
imshow(stackedImage);
width = size(image1, 2);
hold on;
points1=matchedPtsImg1.Location;
points2=matchedPtsImg2.Location;
for i=1:length(theta_stride)
    plot(points1(i, 1), points1(i, 2), 'yo', points2(i, 1) + width, ...
        points2(i, 2), 'y+');
    if theta_stride(i)>0
        line([points1(i, 1) points2(i, 1) + width], [points1(i, 2) points2(i, 2)], ...
            'Color', 'green','linewidth',3);
    else
        line([points1(i, 1) points2(i, 1) + width], [points1(i, 2) points2(i, 2)], ...
            'Color', 'red','linewidth',3);
    end
end
shg