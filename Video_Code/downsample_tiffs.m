% downsample the high-res coronal images and save a movie of downsampled images
% image_path: the folder of high-res coronal images
% save_path: your own folder to save the downsampled tiffs and movie
% scalefactor: scale for downsampling, 0.4 is default, resonable size and good resolution
% FrameRate: frame rate of the movie, default is 5
% intensity_range: intensity limit ranges for rgb channels of the movie, use ImageViewer_LSI to help decide
% put one pair of numbers (0-1) for each channel, i.e. lower and upper bonds
% color_code: numbers in each row are the color code for that channel
% movie_only: if set to 1, will only make a movie without saving downsampled images
%%
image_path = 'V:\Light Sheet Imaging\Suryanarayana_MOS-Br-2_antiGFP-antiRFP-NeuN\raw_coronal_slices';
save_path = 'V:\Light Sheet Imaging\test';
scalefactor = 0.4;
FrameRate = 5;
intensity_range = [0, 0.010695; 0 0.0039; 0 0.056379]; % for channel 1, 2, 3
color_code = [1 1 0;
              1 0 1;
              0 0 1]; % numbers in each row are the color code for that channel
movie_only = 1;

xres = 1.8;
yres = 1.8;
zres = 4;

vidObj = VideoWriter([save_path '\movie.mp4'], 'MPEG-4');
set(vidObj, 'FrameRate', FrameRate, 'Quality', 100);
open(vidObj);
dirlist = dir(image_path);
nimages = 0;
clear image_names imageIDs;
for i = 3:numel(dirlist)
    if numel(dirlist(i).name) > 4
        if strcmp(dirlist(i).name(end-3:end), 'tiff')
            nimages = nimages+1;
            image_names{nimages} = dirlist(i).name;
            imageIDs(nimages) = str2double(dirlist(i).name(end-7:end-5));
        end
    end
end
[~, IDsort] = sort(imageIDs);
image_names = image_names(IDsort);
workbar(0, 'Computing Ongoing...', 'Progress'); 
for i = 1:nimages
    iminfo = imfinfo([image_path '\' image_names{i}], 'tiff');
    for j = 1:numel(iminfo)
        im = imread([image_path '\' image_names{i}], 'tiff', 'Index', j);
        if scalefactor ~= 1
            im = imresize(im, scalefactor);
        end
        im_frame = imresize(im, 'Scale', [1, yres/zres]);
        if j == 1
            if ~movie_only
                imwrite(im, [save_path '\downsampled_' image_names{i}], 'tiff', 'Compression', 'none');
            end
            vidFrame = zeros(size(im_frame, 1), size(im_frame, 2), 3);
        else
            if ~movie_only
                imwrite(im, [save_path '\downsampled_' image_names{i}], 'tiff', 'Compression', 'none', 'WriteMode', 'append');
            end
        end
        im_frame = double(im_frame)/(2^16-1);
        im_frame = imadjust(im_frame, intensity_range(j, :), [0 1]);
        
        for k = 1:3
            vidFrame(:, :, k) = vidFrame(:, :, k)+im_frame*color_code(j, k);
        end
    end
    vidFrame = min(vidFrame, 1);
    vidFrame = max(vidFrame, 0);
    writeVideo(vidObj, vidFrame);
    workbar(i/nimages, [num2str(i) '/' num2str(nimages)], 'Progress'); 
end
close(vidObj);
msgbox('Done !');