close all
clear
clc

% This code is set-up to prompt for input so that it is easier to run.

%% Setting up the code.
% Install the Computer Vision Toolbox on MATLAB in order
% for the set-up to work. 

% Set up video reader -- Video first, and then set the output data
% type to double for higher quality.
prompt_vid = "What is the file name of your video? Don't forget" + ...
    "to use quotation marks." + ...
    "\n";
video = input(prompt_vid);
videoFReader = vision.VideoFileReader(video,'VideoOutputDataType', 'double');

% Set up video player -- This will cause a video player to appear.
vidPlayer    = vision.DeployableVideoPlayer;

% Set up tracker      -- Returns a tracker that tracks an object
% by using the CAMShift algorithm
tracker      = vision.HistogramBasedTracker;

%% Set important information.
prompt_vert = "What is the vertical dimension of your video in " + ...
    "pixels?\n";
vert = input(prompt_vert);

prompt_fps = "What is your camera's FPS (frames per second)?" + ...
    "\n";
fps = input(prompt_fps);

prompt_descend = "Should the image scanner contract downards?" + ...
    " Answer with true or false." + "\n";
descending = input(prompt_descend);

if descending == false
    prompt_solids = "What was the final height of your solids " + ...
        "upon settling?" + "\n";
    solids_height = input(prompt_solids);
else
    solids_height = 0;
end

water_height = 15.5;
vert_ratio = water_height/vert;

% Set up the empty arrays in which height of sediments and the
% confidence of the tracker will be stored.
height = [];
confidence = [];

%% Initialize the tracker.
% Retrieve the first frame from the video and display the image.
img          = step(videoFReader);
figure
imshow(img);

% Select the region of interest.
h            = imrect;
wait(h);
boxLoc       = getLoc(h);
imgYcbcr     = rgb2ycbcr(img);
initializeObject(tracker,imgYcbcr(:,:,2),boxLoc);

% Note: 'Ycbcr' refers to a color space. When we call it, we are
% identifying the color in the chosen area.

%% Track the sediments.
idx = 1;

% Set a loop while the video is playing.
while ~isDone(videoFReader)
    img = step(videoFReader);
    imgYcbcr = rgb2ycbcr(img);
    [sbox,~,score(idx)] = step(tracker,imgYcbcr(:,:,2));

    % sbox refers to the tracker box on the screen. The output is the
    % coordinates of the upper left corner, its width, and its height.

    if score(idx) > 0.5
        out = insertShape(img,'Rectangle',sbox);
        % If confidence is high, display the tracker box and take note of
        % its coordinates.

    else
        % If confidence is low, scan the image again for sediments.  
        boxLoc = segmenting(img,5000);
        if (~isempty(boxLoc))
            % If the sediments are found, reinitialize and start again.
            initializeSearchWindow(tracker,boxLoc)
            [sbox,~,score(idx)] = step(tracker,imgYcbcr(:,:,2));
            out = insertShape(img,'Rectangle',sbox);
        else
            % If the sediment height is not found, output the
            % unaltered frame and take note of the whole frame.
            out = img;
        end
   end
    
    step(vidPlayer,out);       
    idx = idx+1;
    % Display the frame being scanned.

    % Take note of the height every second.
    if mod(idx,fps)==0
        if descending == true
            % If the box is moving normally (with the sediment), take note
            % of its height.

            coord = sbox(4);
        else
            % If the box is moving in the opposite direction, take note of
            % the y-coordinate of the bottom right corner.

            coord = vert - sbox(2) - sbox(4);
        end

    % Convert the height from pixels to centimeters and account for
    % negative readings for height.

        if coord > 0
            height(idx) = coord*vert_ratio;
            confidence(idx) = score(idx-1);
        else
            height(idx) = 0.01;
            confidence(idx) = score(idx-1);
        end
    end
end

%% Cleaning recorded data.
% Clean the height and confidence arrays. The code above only records the
% height/confidence every second -- every other frame is recorded as 0.
% Thus, we remove all 0s from both arrays.

new_height = height(height~=0);
new_confid = confidence(confidence~=0);

% Occasionally, the tracker fluctuates wildly during low confidence
% periods. The following code accounts for these fluctuations so that the
% plot looks more steady.

for i=2:length(new_height)
    if new_height(i) - new_height(i-1) > new_height(i)*(1 - ...
            new_confid(i))*new_confid(i)
        new_height(i) = new_height(i-1);
    else
    end
end

% Close the video players and readers.
release(vidPlayer);
release(videoFReader);

% In cases where descending = false, account for the final solids height.
new_height = new_height + solids_height;

%% Plotting the graph.
xlabel('Time (seconds)')
ylabel('Height (cm)')
figure;
plot(new_height)