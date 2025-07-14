function varargout = ImageViewer_LSI(varargin)
%*---------------------------------------------------------------------*
% Huang Lab
% Cold Spring Harbor Laboratory
% Author : Xu An, June 2017
% xan@cshl.edu
% Version: 1.0
%
% Update 1: Xu An, July 2019
% xan@cshl.edu
% Enable color rendering
%
% Update 2: Xu An, Sep 2019
% xan@cshl.edu
% Enable manual cell marking
%
% Update 3: Xu An, Oct 2021
% xu.an@duke.edu
% Make it compatible with STP data set
%
% Update 4: Xu An, Feb 2024
% xu.an@duke.edu
% Make it compatible with Light-sheet imaging data set
%*---------------------------------------------------------------------*
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @untitled_OpeningFcn, ...
    'gui_OutputFcn',  @untitled_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end

function untitled_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;

handles.radiobutton_selected = 1;
handles.color_code = [1 0 0;
                      0 1 0;
                      0 0 1
                      1 1 0
                      1 0 1
                      0 1 1
                      1 1 1]; % R, G, B, Y, M, C, W
handles.color_code_current = [1 0 0;
                              0 1 0;
                              0 0 1];
set(handles.figure1, 'KeyPressFcn', @quit_marking);
handles.hpoint = [];
guidata(hObject, handles);

function varargout = untitled_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;

function ImageFolder_Callback(hObject, eventdata, handles)
if ~isfield(handles, 'DataDir')
    DataDir = uigetdir('C:\', 'Select an image folder to analyze');
else
    DataDir = uigetdir(handles.DataDir, 'Select an image folder to analyze');
end
if DataDir ~= 0
    handles.DataDir = DataDir;
    set(handles.ImageDir, 'String', DataDir);
    set(handles.image_list, 'String', []);
    Update_Callback(hObject, eventdata, handles);
end

function scalefactor_text_Callback(hObject, eventdata, handles)

function scalefactor_text_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ImageResize_Callback(hObject, eventdata, handles)
curpwd = pwd;
cd(handles.DataDir);
DirList = dir;
n = 0;
for i = 3:size(DirList)
    if strncmp(DirList(i).name(end-3:end), '.jp2', 4)
        n = n+1;
        image_name{n} = DirList(i).name;
    end
end
scalefactor = str2double(get(handles.scalefactor_text, 'String'));
if isnan(scalefactor)
    scalefactor = 0.1;
end
workbar(0, 'Computing Ongoing...', 'Progress'); 
for i = 1:numel(image_name)
    im = imread(image_name{i});
    if scalefactor < 1 && scalefactor > 0
        im = imresize(im, scalefactor);
    end
%     save([image_name{i}(1:end-4) '.mat'], 'im');
    for j = 1:size(im, 3)
        if j == 1
            imwrite(im(:,:,j), [image_name{i}(1:end-4) '.tiff'], 'tiff', 'Compression', 'none');
        else
            imwrite(im(:,:,j), [image_name{i}(1:end-4) '.tiff'], 'tiff', 'Compression', 'none', 'WriteMode', 'append');
        end
    end
    workbar(i/numel(image_name), [num2str(i) '/' num2str(numel(image_name))], 'Progress'); 
end
cd (curpwd);

function Update_Callback(hObject, eventdata, handles)
curpwd = pwd;
cd(handles.DataDir);
DirList = dir;
n = 0;
image_name = [];
for i = 3:size(DirList)
    if strncmp(DirList(i).name(end-4:end), '.tiff', 5)
        n = n+1;
        image_name{n} = DirList(i).name(1:end-5);
        STP = 0;
    end
    
    if strncmp(DirList(i).name(end-3:end), '.jp2', 4) && ~isempty(strfind(DirList(i).name, 'StitchedImage_Z'))
        n = n+1;
        image_name{n} = DirList(i).name(1:end-4);
        STP = 1;
    end
end
if n == 0
    guidata(hObject, handles);
    cd(curpwd);
    return;
end
number_sort = zeros(1, n);
for i = 1:n
    if ~STP
        try
            keywordIDs = strfind(image_name{i}, '_');
            imageid = str2double(image_name{i}(keywordIDs(end)+1:end));
            number_sort(i) = imageid;
        catch
            number_sort(i) = NaN;
        end
    else
        number_sort(i) = str2double(image_name{i}(16:18));
    end
end
[~, index] = sort(number_sort);
image_name = image_name(index);
set(handles.image_list, 'String', image_name);
set(handles.image_list, 'Value', []);
set(handles.colorchannel_list, 'String', []);
set(handles.colorchannel_list, 'Value', []);
try
    load('Adjustment_Parameters.mat');
    if ~exist('area_ratios', 'var')
        area_ratios = ones(1, n)*0.2;
    end
    if ~exist('thrs', 'var')
        thrs = zeros(n, 3)/0;
    end
    if ~exist('cell_detection', 'var')
        cell_detection = cell(3, n);
    end
    save('Adjustment_Parameters.mat', 'rotation_angles', 'mirror_image', 'display_ranges', 'area_ratios', 'thrs', 'cell_detection');
catch
    rotation_angles = zeros(1, n);
    mirror_image = zeros(1, n);
    display_ranges = zeros(3, 2, n);
    display_ranges(:, 2, :) = 1;
    area_ratios = ones(1, n)*0.2;
    thrs = zeros(n, 3)/0;
    cell_detection = cell(3, n);
    save('Adjustment_Parameters.mat', 'rotation_angles', 'mirror_image', 'display_ranges', 'area_ratios', 'thrs', 'cell_detection');
end
handles.rotation_angles = rotation_angles;
handles.mirror_image = mirror_image;
handles.display_ranges = display_ranges;
handles.area_ratios = area_ratios;
handles.thrs = thrs;
handles.cell_detection = cell_detection;
handles.STP = STP;
guidata(hObject, handles);
cd(curpwd);

function image_list_Callback(hObject, eventdata, handles)
if numel(get(handles.image_list, 'Value')) ~= 1
    return;
end
image_name = get(handles.image_list, 'String');
image_name = image_name{get(handles.image_list, 'Value')};
curpwd = pwd;
cd(handles.DataDir);
if get(handles.highres_check, 'Value')
    im = imread([image_name '.jp2']);
    msgbox('Image has been loaded!', '')
else
    try
        iminfo = imfinfo([image_name '.tiff'], 'tiff');
        im = zeros(iminfo(1).Height, iminfo(1).Width, numel(iminfo), 'uint16');
        for i = 1:numel(iminfo)
            im(:,:,i) = imread([image_name '.tiff'], 'tiff', 'Index', i);
        end
    catch
        load([image_name '.mat']);
    end
end
handles.im = im;
handles.Restrict_button.Value = 0;
guidata(hObject, handles);
channel = get(handles.colorchannel_list, 'Value');
if isempty(channel)
    set(handles.colorchannel_list, 'Value', 1);
end
color_channel = cell(1, size(im, 3));
for i = 1:size(im, 3)
    color_channel{i} = ['Channel ' num2str(i)];
end
if max(channel) > size(im, 3)
    set(handles.colorchannel_list, 'Value', []);
    set(handles.colorchannel_list, 'String', color_channel);
    cd(curpwd);
    return;
end
set(handles.colorchannel_list, 'String', color_channel);
colorchannel_list_Callback(hObject, eventdata, handles);
intensity_histogram(hObject, eventdata, handles);
rotation_angle = handles.rotation_angles(get(handles.image_list, 'Value'));
set(handles.rotation_text, 'String', num2str(rotation_angle));
set(handles.MirrorImage_check, 'Value', handles.mirror_image(get(handles.image_list, 'Value')));
area_ratio = handles.area_ratios(get(handles.image_list, 'Value'));
set(handles.arearatio_text, 'String', num2str(area_ratio));
cd(curpwd);

function intensity_histogram(hObject, eventdata, handles)
if handles.Restrict_button.Value
    return;
end
imageid = get(handles.image_list, 'Value');
if numel(imageid) ~= 1
    errordlg('Please choose only 1 image!');
    return;
end
channel = get(handles.colorchannel_list, 'Value');
cla(handles.imhist_axes);
im = handles.im;
if get(handles.radiobutton1, 'Value')
    if ~isempty(intersect(1, channel))
        im = im(:, :, 1);
        channel = 1;
    else
        set(handles.imhist_axes, 'UserData', []);
        return;
    end
end
if  get(handles.radiobutton2, 'Value')
    if ~isempty(intersect(2, channel))
        im = im(:, :, 2);
        channel = 2;
    else
        set(handles.imhist_axes, 'UserData', []);
        return;
    end
end
if  get(handles.radiobutton3, 'Value')
    if ~isempty(intersect(3, channel))
        im = im(:, :, 3);
        channel = 3;
    else
        set(handles.imhist_axes, 'UserData', []);
        return;
    end
end
hcolor = handles.color_code_current(channel, :);
if ~any(hcolor ~= 1)
    hcolor = hcolor*0.8;
end
im = double(im)/(2^16-1);
if handles.rotation_angles(imageid) ~= 0
    im = imrotate(im, handles.rotation_angles(imageid), 'bicubic');
end
if handles.mirror_image(imageid) ~= 0
    im = fliplr(im);
end
[count, x] = imhist(im, 2^16);
x = x(min(find(count ~= 0)):max(find(count ~= 0)));
count = count(min(find(count ~= 0)):max(find(count ~= 0)));
stem(handles.imhist_axes, x, count, 'Color', hcolor, 'Marker', 'none');
try
    ylim(handles.imhist_axes, [0 max(count)]);
    xlim(handles.imhist_axes, [x(1) x(end)]);
end
hold(handles.imhist_axes, 'on');
hlow = plot(handles.imhist_axes, [max(handles.display_ranges(channel, 1, imageid), x(1)) max(handles.display_ranges(channel, 1, imageid), x(1))], [0 max(count)], '--', 'Color', [0 0 0], 'LineWidth', 1);
hhigh = plot(handles.imhist_axes, [min(handles.display_ranges(channel, 2, imageid), x(end)) min(handles.display_ranges(channel, 2, imageid), x(end))], [0 max(count)], '--', 'Color', [0 0 0], 'LineWidth', 1);
set(handles.imhist_axes, 'UserData', [hlow hhigh]);
title(handles.imhist_axes, ['Image ' num2str(imageid)]);
set(handles.MinValue_slider, 'Min', x(1));
set(handles.MinValue_slider, 'Max', x(end));
set(handles.MinValue_slider, 'SliderStep', [1/2^8 1/2^8]);
set(handles.MinValue_slider, 'Value', max(handles.display_ranges(channel, 1, imageid), x(1)));
set(handles.minvalue_text, 'String', num2str(max(handles.display_ranges(channel, 1, imageid), x(1))));
set(handles.MaxValue_slider, 'Min', x(1));
set(handles.MaxValue_slider, 'Max', x(end));
set(handles.MaxValue_slider, 'SliderStep', [1/2^8 1/2^8]);
set(handles.MaxValue_slider, 'Value', min(handles.display_ranges(channel, 2, imageid), x(end)));
set(handles.maxvalue_text, 'String', num2str(min(handles.display_ranges(channel, 2, imageid), x(end))));
clear im;
    
function image_list_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function colorchannel_list_Callback(hObject, eventdata, handles)
persistent hrect
try
    delete(hrect);
end
imageid = get(handles.image_list, 'Value');
if numel(imageid) ~= 1
    errordlg('Please choose only 1 image!');
    return;
end
channel = get(handles.colorchannel_list, 'Value');
im = handles.im;
for i = 1:numel(channel)
    temp = im(:, :, channel(i));
    if channel(i) ~= 4
        temp = double(temp)/(2^16-1);
        if handles.rotation_angles(imageid) ~= 0
            temp = imrotate(temp, handles.rotation_angles(imageid), 'bicubic');
        end
        if handles.mirror_image(imageid) ~= 0
            temp = fliplr(temp);
        end
        if get(handles.BBox_check, 'Value')
            mask = temp > handles.display_ranges(channel(i), 2, imageid);
        end
        temp = imadjust(temp, [handles.display_ranges(channel(i), 1, imageid) handles.display_ranges(channel(i), 2, imageid)], [0 1]);
        if i == 1
            im_display = zeros([size(temp), 3]);
        end
        for j = 1:3
            im_display(:, :, j) = im_display(:, :, j)+temp*handles.color_code_current(channel(i), j);
        end
        eval(['set(handles.radiobutton' num2str(channel(i)) ', ''Visible'', ''on'');']);
    else
        if handles.rotation_angles(imageid) ~= 0
            temp = imrotate(temp, handles.rotation_angles(imageid), 'nearest');
        end
        if handles.mirror_image(imageid) ~= 0
            temp = fliplr(temp);
        end
        if i == 1
            im_display = zeros([size(temp), 3]);
        end
        for j = 1:3
            temp_im_display = im_display(:, :, j);
            temp_im_display(temp == 255) = 1;
            im_display(:, :, j) = temp_im_display;
        end
    end
end
im_display = min(im_display, 1);
try
    handles.image_axes.Children.CData = im_display;
catch
    cla(handles.image_axes);
    image(im_display, 'Parent', handles.image_axes);
    axis(handles.image_axes, 'image');
    axis(handles.image_axes, 'off');
end
if get(handles.fix_distortion_check, 'Value')
    set(handles.image_axes, 'DataAspectRatio', [4, 1.8, 1.8]);
else
    set(handles.image_axes, 'DataAspectRatio', [1, 1, 1]);
end
if numel(channel) == 1
    if channel == 1
        set(handles.radiobutton1, 'Value', 1);
    end
    if channel == 2
        set(handles.radiobutton2, 'Value', 1);
    end
    if channel == 3
        set(handles.radiobutton3, 'Value', 1);
    end
    if ~strcmp(hObject.Style, 'slider')
        handles.Restrict_button.Value = 0;
        intensity_histogram(hObject, eventdata, handles);
    end
    
    if get(handles.BBox_check, 'Value')
        se = strel('disk', 3);
        mask = imerode(mask, se);
        mask = imdilate(mask, se);
        stats = regionprops(mask, 'Area', 'BoundingBox');
        if numel(stats) ~= 0
            area = zeros(1, numel(stats));
            boundingbox = zeros(numel(stats), 4);
            for j = 1:numel(stats)
                area(j) = stats(j).Area;
                boundingbox(j, :) = stats(j).BoundingBox;
                boundingbox(j, 3:4) = boundingbox(j, 1:2)+boundingbox(j, 3:4);
            end
            maxarea = max(area);
            boundingbox(area < maxarea*handles.area_ratios(imageid), :) = [];
            if size(boundingbox, 1) ~= 1
                boundingbox_all(1:2) = min(boundingbox(:, 1:2));
                boundingbox_all(3:4) = max(boundingbox(:, 3:4));
            else
                boundingbox_all = boundingbox;
            end
            boundingbox_all(3:4) = boundingbox_all(3:4)-boundingbox_all(1:2);
            hold(handles.image_axes, 'on');
            hrect = rectangle('Position', boundingbox_all, 'EdgeColor', [1 1 1], 'LineWidth', 1, 'Parent', handles.image_axes);
        end
    end
end
off_channel = setdiff([1 2 3], channel);
if ~isempty(off_channel)
    for i = 1:numel(off_channel)
        eval(['set(handles.radiobutton' num2str(off_channel(i)) ', ''Visible'', ''off'');']);
    end
end
clear im temp im_display;

function colorchannel_list_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function MaxValue_slider_Callback(hObject, eventdata, handles)
imageid = get(handles.image_list, 'Value');
if numel(imageid) ~= 1
    errordlg('Please choose only 1 image!');
    return;
end
maxvalue = get(handles.MaxValue_slider, 'Value');
minvalue = get(handles.MinValue_slider, 'Value');
if maxvalue <= minvalue
    maxvalue = minvalue+1/2^8;
    set(handles.MaxValue_slider, 'Value', maxvalue);
end
set(handles.maxvalue_text, 'String', num2str(maxvalue));
temp = get(handles.imhist_axes, 'UserData');
if isempty(temp)
    return;
end
hhigh = temp(2);
set(hhigh, 'XData', [maxvalue maxvalue]);
if get(handles.radiobutton1, 'Value')
    channel = 1;
end
if get(handles.radiobutton2, 'Value')
    channel = 2;
end
if get(handles.radiobutton3, 'Value')
    channel = 3;
end
handles.display_ranges(channel, 2, imageid) = maxvalue;
update_parameters(hObject, eventdata, handles);
colorchannel_list_Callback(hObject, eventdata, handles);

function MaxValue_slider_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function MinValue_slider_Callback(hObject, eventdata, handles)
imageid = get(handles.image_list, 'Value');
if numel(imageid) ~= 1
    errordlg('Please choose only 1 image!');
    return;
end
maxvalue = get(handles.MaxValue_slider, 'Value');
minvalue = get(handles.MinValue_slider, 'Value');
if minvalue >= maxvalue
    minvalue = maxvalue-1/2^8;
    set(handles.MinValue_slider, 'Value', minvalue);
end
set(handles.minvalue_text, 'String', num2str(minvalue));
temp = get(handles.imhist_axes, 'UserData');
if isempty(temp)
    return;
end
hlow = temp(1);
set(hlow, 'XData', [minvalue minvalue]);
if get(handles.radiobutton1, 'Value')
    channel = 1;
end
if get(handles.radiobutton2, 'Value')
    channel = 2;
end
if get(handles.radiobutton3, 'Value')
    channel = 3;
end
if ~get(handles.BBox_check, 'Value')
    handles.display_ranges(channel, 1, imageid) = minvalue;
    update_parameters(hObject, eventdata, handles);
end
colorchannel_list_Callback(hObject, eventdata, handles);

function MinValue_slider_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function rotation_text_Callback(hObject, eventdata, handles)
imageid = get(handles.image_list, 'Value');
if numel(imageid) ~= 1
    errordlg('Please choose only 1 image!');
    return;
end
rotation_angle = str2double(get(handles.rotation_text, 'String'));
handles.rotation_angles(imageid) = rotation_angle;
set(handles.BBox_check, 'Value', 0);
update_parameters(hObject, eventdata, handles);
colorchannel_list_Callback(hObject, eventdata, handles);
im_display = handles.image_axes.Children.CData;
[height, width, ~] = size(im_display);
im_display(round(height/2-height/500/2):round(height/2+height/500/2), :, :) = 1;
im_display(:, round(width/2-width/500/2):round(width/2+width/500/2), :) = 1;
handles.image_axes.Children.CData = im_display;
intensity_histogram(hObject, eventdata, handles);

function rotation_text_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function RotateImage_Callback(hObject, eventdata, handles)

function ImageDir_Callback(hObject, eventdata, handles)

function ImageDir_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function radiobutton1_Callback(hObject, eventdata, handles)
handles.radiobutton_selected = 1;
handles.Restrict_button.Value = 0;
guidata(hObject, handles);
intensity_histogram(hObject, eventdata, handles);

function radiobutton2_Callback(hObject, eventdata, handles)
handles.radiobutton_selected = 2;
handles.Restrict_button.Value = 0;
guidata(hObject, handles);
intensity_histogram(hObject, eventdata, handles);

function radiobutton3_Callback(hObject, eventdata, handles)
handles.radiobutton_selected = 3;
handles.Restrict_button.Value = 0;
guidata(hObject, handles);
intensity_histogram(hObject, eventdata, handles);

function maxvalue_text_Callback(hObject, eventdata, handles)
maxvalue = str2double(get(handles.maxvalue_text, 'String'));
if maxvalue > get(handles.MaxValue_slider, 'Max')
    maxvalue = get(handles.MaxValue_slider, 'Max');
    set(handles.maxvalue_text, 'String', num2str(maxvalue));
end
if maxvalue < get(handles.MaxValue_slider, 'Min')
    maxvalue = get(handles.MaxValue_slider, 'Min');
    set(handles.maxvalue_text, 'String', num2str(maxvalue));
end
set(handles.MaxValue_slider, 'Value', maxvalue);
MaxValue_slider_Callback(hObject, eventdata, handles);

function maxvalue_text_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function minvalue_text_Callback(hObject, eventdata, handles)
minvalue = str2double(get(handles.minvalue_text, 'String'));
if minvalue < get(handles.MinValue_slider, 'Min')
    minvalue = get(handles.MinValue_slider, 'Min');
    set(handles.minvalue_text, 'String', num2str(minvalue));
end
if minvalue > get(handles.MinValue_slider, 'Max')
    minvalue = get(handles.MinValue_slider, 'Max');
    set(handles.minvalue_text, 'String', num2str(minvalue));
end
set(handles.MinValue_slider, 'Value', minvalue);
MinValue_slider_Callback(hObject, eventdata, handles);

function minvalue_text_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function highres_check_Callback(hObject, eventdata, handles)
image_list_Callback(hObject, eventdata, handles);

function SelectROI_Callback(hObject, eventdata, handles)
if numel(get(handles.image_list, 'Value')) ~= 1
    errordlg('Please choose only 1 image!');
    return;
end
himrect = imrect(handles.image_axes);
fcn = makeConstrainToRectFcn('imrect', get(handles.image_axes,'XLim'),get(handles.image_axes, 'YLim'));
setPositionConstraintFcn(himrect, fcn);
position = wait(himrect);
position = round(position);
delete(himrect);
hold(handles.image_axes, 'on');
rectangle('Position', position, 'Parent', handles.image_axes, 'EdgeColor', [1 1 1], 'LineWidth', 1);
handles.position = position;
ShowCropImage(handles);

function figure1_CreateFcn(hObject, eventdata, handles)

function SaveImage_Callback(hObject, eventdata, handles)
fcolor = get(handles.figure1, 'Color');
set(handles.figure1, 'Color', [0 0 0]);
im = getframe(handles.image_axes);
set(handles.figure1, 'Color', fcolor);
curpwd = pwd;
try
    cd(handles.masterhandles.DataDir);
end
[file, path] = uiputfile('*.tiff','Save');
if file ~= 0
    imwrite(im.cdata, [path file], 'TIFF', 'Compression', 'none');
end
if handles.STP
    if numel(get(handles.image_list, 'Value')) ~= 1
        return;
    end
    image_name = get(handles.image_list, 'String');
    image_name = image_name{get(handles.image_list, 'Value')};
    if ~isfile([handles.DataDir '\' image_name '.tiff'])
        imwrite(handles.im, [handles.DataDir '\' image_name '.tiff'], 'TIFF', 'Compression', 'none');
        msgbox('Done !');
    end
end
cd(curpwd);

function Magic_Callback(hObject, eventdata, handles)
imageIDs = get(handles.image_list, 'Value');
image_names = get(handles.image_list, 'String');
channel = get(handles.colorchannel_list, 'Value');
scalefactor = str2double(get(handles.scalefactor_text, 'String'));
curpwd = pwd;
cd(handles.DataDir);

[file, path] = uiputfile('*.mp4','Save');
if file == 0
    cd(curpwd);
    return;
end

FrameRate = 5;
vidObj = VideoWriter([path file], 'MPEG-4');
set(vidObj, 'FrameRate', FrameRate, 'Quality', 100);
open(vidObj);

workbar(0, 'Computing Ongoing...', 'Progress'); 
for i = 1:numel(imageIDs)
    image_name = image_names{imageIDs(i)};
    disp(image_name);
    if i == 1
        iminfo = imfinfo([image_name '.tiff'], 'tiff');
    end
    im = zeros(iminfo(1).Height, iminfo(1).Width, numel(iminfo), 'uint16');
    for j = 1:numel(iminfo)
        im(:, :, j) = imread([image_name '.tiff'], 'tiff', 'Index', j);
    end

    for j = 1:numel(channel)
        temp = im(:, :, channel(j));
        if channel(j) ~= 4
            temp = double(temp)/(2^16-1);
            temp = imadjust(temp, [handles.display_ranges(channel(j), 1, imageIDs(i)) handles.display_ranges(channel(j), 2, imageIDs(i))], [0 1]);
            if j == 1
                im_display = zeros([size(temp), 3]);
            end
            for k = 1:3
                im_display(:, :, k) = im_display(:, :, k)+temp*handles.color_code_current(channel(j), k);
            end
        else
            if j == 1
                im_display = zeros([size(temp), 3]);
            end
            for k = 1:3
                temp_im_display = im_display(:, :, k);
                temp_im_display(temp == 255) = 1;
                im_display(:, :, k) = temp_im_display;
            end
        end
    end
    if get(handles.fix_distortion_check, 'Value')
        im_display = imresize(im_display, round([iminfo(1).Height*4/1.8 iminfo(1).Width]*scalefactor), 'bicubic');
    else
        if scalefactor ~= 1
            im_display = imresize(im_display, round([iminfo(1).Height iminfo(1).Width]*scalefactor), 'bicubic');
        end
    end
    im_display(im_display < 0) = 0;
    im_display(im_display > 1) = 1;

    writeVideo(vidObj, im_display);
    workbar(i/numel(imageIDs), 'Computing Ongoing...', 'Progress'); 
end
clear im temp im_display;
cd(curpwd);
close(vidObj);
disp('Done!');

function MirrorImage_check_Callback(hObject, eventdata, handles)
imageid = get(handles.image_list, 'Value');
if numel(imageid) ~= 1
    errordlg('Please choose only 1 image!');
    return;
end
handles.mirror_image(imageid) = get(handles.MirrorImage_check, 'Value');
set(handles.BBox_check, 'Value', 0);
update_parameters(hObject, eventdata, handles);
colorchannel_list_Callback(hObject, eventdata, handles);

function arearatio_text_Callback(hObject, eventdata, handles)
imageid = get(handles.image_list, 'Value');
if numel(imageid) ~= 1
    errordlg('Please choose only 1 image!');
    return;
end
area_ratio = str2double(get(handles.arearatio_text, 'String'));
handles.area_ratios(imageid) = area_ratio;
set(handles.BBox_check, 'Value', 1);
update_parameters(hObject, eventdata, handles);
colorchannel_list_Callback(hObject, eventdata, handles);

function arearatio_text_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function BBox_check_Callback(hObject, eventdata, handles)
imageid = get(handles.image_list, 'Value');
if numel(imageid) ~= 1
    set(handles.BBox_check, 'Value', 0);
    errordlg('Please choose only 1 image!');
    return;
end
channel = get(handles.colorchannel_list, 'Value');
if numel(channel) ~= 1
    set(handles.BBox_check, 'Value', 0);
    errordlg('Please choose only 1 Channel!');
    return;
end
if get(handles.BBox_check, 'Value')
    maxvalue = get(handles.MaxValue_slider, 'Value');
    set(handles.BBox_check, 'UserData', maxvalue);
    set(handles.image_list, 'Enable', 'off');
    set(handles.colorchannel_list, 'Enable', 'off');
    set(handles.MinValue_slider, 'Enable', 'off');
    set(handles.minvalue_text, 'Enable', 'off');
    colorchannel_list_Callback(hObject, eventdata, handles);
else
    set(handles.image_list, 'Enable', 'on');
    set(handles.colorchannel_list, 'Enable', 'on');
    set(handles.MinValue_slider, 'Enable', 'on');
    set(handles.minvalue_text, 'Enable', 'on');
    maxvalue = get(handles.BBox_check, 'UserData');
    set(handles.MaxValue_slider, 'Value', maxvalue);
    MaxValue_slider_Callback(hObject, eventdata, handles);
end

function WriteThr_Callback(hObject, eventdata, handles)
imageid = get(handles.image_list, 'Value');
if numel(imageid) ~= 1
    errordlg('Please choose only 1 image!');
    return;
end
channel = get(handles.colorchannel_list, 'Value');
if numel(channel) ~= 1
    errordlg('Please choose only 1 channel!');
    return;
end
maxvalue = str2double(get(handles.maxvalue_text, 'String'));
handles.thrs(imageid, channel) = maxvalue;
update_parameters(hObject, eventdata, handles);
msgbox('Done!');

function update_parameters(hObject, eventdata, handles)
guidata(hObject, handles);
display_ranges = handles.display_ranges;
rotation_angles = handles.rotation_angles;
mirror_image = handles.mirror_image;
area_ratios = handles.area_ratios;
thrs = handles.thrs;
cell_detection = handles.cell_detection;
curpwd = pwd;
cd(handles.DataDir);
save('Adjustment_Parameters.mat', 'rotation_angles', 'display_ranges', 'mirror_image', 'area_ratios', 'thrs', 'cell_detection');
cd(curpwd);

function Range_Apply_All_button_Callback(hObject, eventdata, handles)
maxvalue = str2double(get(handles.maxvalue_text, 'String'));
minvalue = str2double(get(handles.minvalue_text, 'String'));
if get(handles.radiobutton1, 'Value')
    channel = 1;
end
if get(handles.radiobutton2, 'Value')
    channel = 2;
end
if get(handles.radiobutton3, 'Value')
    channel = 3;
end
handles.display_ranges(channel, 1, :) = minvalue;
handles.display_ranges(channel, 2, :) = maxvalue;
update_parameters(hObject, eventdata, handles);

function Area_Apply_All_button_Callback(hObject, eventdata, handles)
area_ratios = str2double(get(handles.arearatio_text, 'String'));
handles.area_ratios(:) = area_ratios;
update_parameters(hObject, eventdata, handles);

function popupmenu1_Callback(hObject, eventdata, handles)
handles.color_code_current(1, :) = handles.color_code(get(handles.popupmenu1, 'Value'), :);
guidata(hObject, handles);
colorchannel_list_Callback(hObject, eventdata, handles);
intensity_histogram(hObject, eventdata, handles);

function popupmenu1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popupmenu2_Callback(hObject, eventdata, handles)
handles.color_code_current(2, :) = handles.color_code(get(handles.popupmenu2, 'Value'), :);
guidata(hObject, handles);
colorchannel_list_Callback(hObject, eventdata, handles);
intensity_histogram(hObject, eventdata, handles);

function popupmenu2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popupmenu3_Callback(hObject, eventdata, handles)
handles.color_code_current(3, :) = handles.color_code(get(handles.popupmenu3, 'Value'), :);
guidata(hObject, handles);
colorchannel_list_Callback(hObject, eventdata, handles);
intensity_histogram(hObject, eventdata, handles);

function popupmenu4_Callback(hObject, eventdata, handles)
handles.color_code_current(4, :) = handles.color_code(get(handles.popupmenu4, 'Value'), :);
guidata(hObject, handles);
colorchannel_list_Callback(hObject, eventdata, handles);
intensity_histogram(hObject, eventdata, handles);

function popupmenu4_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popupmenu3_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Mark_Cells_Callback(hObject, eventdata, handles)
channel = get(handles.colorchannel_list, 'Value');
if numel(channel) ~= 1
    errordlg('Please choose only 1 Channel!');
    return;
end
handles.figure1.UserData = 1;
while handles.figure1.UserData
    h = drawpoint(handles.image_axes, 'Color', [1 1 1]);
    handles.hpoint = [handles.hpoint h];
end
guidata(hObject, handles);

function quit_marking(hObject, eventdata, handles)
if isequal(eventdata.Key, 'escape')
    hObject.UserData = 0;
end

function Save_Positions_Callback(hObject, eventdata, handles)
channel = get(handles.colorchannel_list, 'Value');
if numel(channel) ~= 1
    errordlg('Please choose only 1 Channel!');
    return;
end
imageid = get(handles.image_list, 'Value');
if numel(imageid) ~= 1
    errordlg('Please choose only 1 image!');
    return;
end
temp = [];
for i = 1:numel(handles.hpoint)
    if isvalid(handles.hpoint(i))
        temp = [temp; handles.hpoint(i).Position];
    end
end
handles.cell_detection{channel, imageid} = temp;
update_parameters(hObject, eventdata, handles);
delete(handles.hpoint);
handles.hpoint = [];
guidata(hObject, handles);
msgbox([num2str(size(temp, 1)) ' Cell Position Saved']);

function Recover_Markers_Callback(hObject, eventdata, handles)
channel = get(handles.colorchannel_list, 'Value');
if numel(channel) ~= 1
    errordlg('Please choose only 1 Channel!');
    return;
end
imageid = get(handles.image_list, 'Value');
if numel(imageid) ~= 1
    errordlg('Please choose only 1 image!');
    return;
end
if ~isempty(handles.cell_detection{channel, imageid})
    handles.hpoint = images.roi.Point.empty;
    temp = handles.cell_detection{channel, imageid};
    for i = 1:size(temp, 1)
        handles.hpoint(i) = drawpoint(handles.image_axes, 'Color', [1 1 1], 'Position', temp(i, :));
    end
end
guidata(hObject, handles);

function fix_distortion_check_Callback(hObject, eventdata, handles)
if get(handles.fix_distortion_check, 'Value')
    set(handles.image_axes, 'DataAspectRatio', [4, 1.8, 1.8]);
else
    set(handles.image_axes, 'DataAspectRatio', [1, 1, 1]);
end

function Restrict_button_Callback(hObject, eventdata, handles)
if handles.Restrict_button.Value
    maxvalue = get(handles.MaxValue_slider, 'Value');
    minvalue = get(handles.MinValue_slider, 'Value');
    set(handles.MinValue_slider, 'Min', minvalue);
    set(handles.MinValue_slider, 'Max', maxvalue);
    set(handles.MaxValue_slider, 'Min', minvalue);
    set(handles.MaxValue_slider, 'Max', maxvalue);
    set(handles.imhist_axes, 'XLim', [minvalue maxvalue]);
else
    intensity_histogram(hObject, eventdata, handles);
end
