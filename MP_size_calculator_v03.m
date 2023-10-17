%% This is used to estimate size of MP particles
% Inputs are black-white png images with particles in white color

% 2023/07/18   By Zijiang Yang
% 2023/07/18b  Add more parameters for image
% 2023/07/19   Fix MaxFeretDiameter bugs


%% Prepare workspace
clc, clear, close all



%% Input parameters
tic
% Image filename
filename = 'Box 01.TIF';

% filename2 = 'Slide22.TIF.tif';
% image2 = imread(filename2);

% Factor for conversion from pixels to physical units
f = 176.88; % replace with your actual factor

% Filter out objects smaller than x pixels: shreshold
sh = 50;


%% Image processing
% Load Image
image = imread(filename);



% Convert to grayscale if image is color
if ndims(image) == 3
    image = rgb2gray(image);
end

% Convert the image to binary
binaryImage = imbinarize(image);

% Identify the connected components in the image
connectedComponents = bwconncomp(binaryImage, 8); % Using 8-connectivity

% Measure the area, major axis length, minor axis length, and perimeter of each object in pixels
stats = regionprops(connectedComponents, 'Area', 'MajorAxisLength', 'MinorAxisLength', 'Centroid', 'Orientation', 'Perimeter','MaxFeretProperties');

% Filter out objects smaller than 10 pixels
stats([stats.Area] < sh) = [];

stats2 = stats;
% Remove un-wanted MaxFeretProperties
newStats2 = struct('Area', [], 'Centroid', [], 'MajorAxisLength', [], 'MinorAxisLength', [], 'Orientation', [], 'Perimeter', [], 'MaxFeretDiameter', []);

for i = 1:numel(stats2)
    fieldsToRemove = {'MaxFeretAngle', 'MaxFeretCoordinates'};
    fieldsToCopy = setdiff(fieldnames(stats2(i)), fieldsToRemove);
    for j = 1:numel(fieldsToCopy)
        newStats2(i).(fieldsToCopy{j}) = stats2(i).(fieldsToCopy{j});
    end
end

stats = newStats2;


%% Calculate the additional properties and convert the measurements to physical units
for idx = 1:length(stats)

    stats(idx).Index = idx;
    stats(idx).Area = stats(idx).Area / (f^2);
    stats(idx).MajorAxisLength = stats(idx).MajorAxisLength / f;
    stats(idx).MinorAxisLength = stats(idx).MinorAxisLength / f;
    stats(idx).Perimeter = stats(idx).Perimeter / f;
    
    % Calculating Circularity (Roundness or Compactness)
    stats(idx).Circularity = 4*pi*stats(idx).Area / (stats(idx).Perimeter^2);
    
    
    % Calculate max Feret diameter
    stats(idx).MaxFeretDiameter = stats(idx).MaxFeretDiameter / f;
end

% Create a table to store the measurements
objectTable = struct2table(stats);

% % Display the original image
% figure;
% imshow(image);

%% Display the original image
figure;
imshow(binaryImage);

% Changes the colormap of the image to display white background and light gray objects
colormap([1 1 1; 0.9 0.9 0.9]);

hold on;

% Loop through the stats array and overlay major and minor axis on the image
for idx = 1:length(stats)
    % Overlay the major and minor axis lengths
    str = sprintf('#%d, (%.2f, %.2f, %.2f)', stats(idx).Index, stats(idx).MajorAxisLength, stats(idx).MinorAxisLength,stats(idx).MaxFeretDiameter); % Index added in the string <------ HERE is the addition
    text(stats(idx).Centroid(1), stats(idx).Centroid(2), str, 'Color', 'r', 'FontSize', 10);
end


% Set the title of the image
[~, name, ~] = fileparts(filename);
title(sprintf('%s, Number of Objects: %d, (lmj, lmn, lmax) unit: mm', name, numel(stats)));


%% Display the original image (index only)
figure;
imshow(binaryImage);

% Changes the colormap of the image to display white background and light gray objects
colormap([1 1 1;  0.9 0.9 0.9]);

hold on;

% Loop through the stats array and overlay major and minor axis on the image
for idx = 1:length(stats)
    % Overlay the major and minor axis lengths
    str = sprintf('#%d', stats(idx).Index); % Index added in the string <------ HERE is the addition
    text(stats(idx).Centroid(1), stats(idx).Centroid(2), str, 'Color', 'r', 'FontSize', 10);
end


% Set the title of the image
[~, name, ~] = fileparts(filename);
title(sprintf('%s, Number of Objects: %d', name, numel(stats)));

%disp(objectTable)
% Add the .csv extension to the filename
filename = strcat(filename, '.csv');

% Export the table as a CSV file
writetable(objectTable, filename);

toc



