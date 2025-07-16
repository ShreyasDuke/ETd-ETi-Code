% resection of raw horizontal sections in coronal orientation
% keep raw resolution during resection, will take a long time to run
% data_path: home folder for the sample (the folder for all channels), e.g. G:\LifeCanvas CRO datasets\Suryanarayana_V5_antiGFP-antiRFP-NeuN
% save_path: your own folder to save the high-res tiffs
% interval: thickness of resection, default is 54 microns
%%
data_path = {'G:\⁭LifeCanvas CRO datasets\RawData_ShreyasSuryanarayana_012225_Duke\Suryanarayana_PONS-Br-4_antiGFP-antiRFP-NeuN',...
    'G:\⁭LifeCanvas CRO datasets\RawData_ShreyasSuryanarayana_012225_Duke\Suryanarayana_RetOrb-Br-1_antiGFP-antiRFP-NeuN',...
    'G:\⁭LifeCanvas CRO datasets\RawData_ShreyasSuryanarayana_012225_Duke\Suryanarayana_RetOrb-Br-2_antiGFP-antiRFP-NeuN',...
    'G:\⁭LifeCanvas CRO datasets\RawData_ShreyasSuryanarayana_012225_Duke\Suryanarayana_SpCor-Br-1_antiGFP-antiRFP-NeuN',...
    'G:\⁭LifeCanvas CRO datasets\RawData_ShreyasSuryanarayana_012225_Duke\Suryanarayana_SpCor-Br-2_antiGFP-antiRFP-NeuN'};
save_path = {'V:\Light Sheet Imaging\Suryanarayana_PONS-Br-4_antiGFP-antiRFP-NeuN\raw_coronal_slices',...
    'V:\Light Sheet Imaging\Suryanarayana_RetOrb-Br-1_antiGFP-antiRFP-NeuN\raw_coronal_slices',...
    'V:\Light Sheet Imaging\Suryanarayana_RetOrb-Br-2_antiGFP-antiRFP-NeuN\raw_coronal_slices',...
    'V:\Light Sheet Imaging\Suryanarayana_SpCor-Br-1_antiGFP-antiRFP-NeuN\raw_coronal_slices',...
    'V:\Light Sheet Imaging\Suryanarayana_SpCor-Br-2_antiGFP-antiRFP-NeuN\raw_coronal_slices'};
xres = 1.8;
yres = 1.8;
zres = 4;
interval = 54; % microns default 54
interval_pixels = interval/xres;

workbar(0, 'Computing Ongoing...', 'Progress');
for i = 1:numel(data_path)
    dirlist = dir(data_path{i});
    nchannels = 0;
    clear channel_names section_numbers;
    for j = 3:size(dirlist)
        if dirlist(j).isdir
            nchannels = nchannels+1;
            channel_names{nchannels} = dirlist(j).name;
        end
    end

    for j = 1:nchannels
        channel_path = [data_path{i} filesep channel_names{j}];
        dirlist = dir(channel_path);
        nhorizontal = 0;
        image_names = cell(0);
        for k = 3:numel(dirlist)
            if numel(dirlist(k).name) > 4
                if strcmp(dirlist(k).name(end-2:end), 'tif')
                    nhorizontal = nhorizontal+1;
                    image_names{nhorizontal} = dirlist(k).name;
                    ids = strfind(image_names{nhorizontal}, '_');
                    if numel(ids) == 3
                        section_numbers(nhorizontal) = str2num(image_names{nhorizontal}(ids(2)+1:ids(3)-1));
                    elseif numel(ids) == 2
                        section_numbers(nhorizontal) = str2num(image_names{nhorizontal}(ids(2)+1:end-4));
                    end
                end
            end
        end
        [~, IDsort] = sort(section_numbers);
        image_names = image_names(IDsort);

        iminfo = imfinfo([channel_path filesep image_names{1}], 'tiff');
        nrow = nhorizontal;
        ncolumn = iminfo.Width;
        ncoronal = ceil(iminfo.Height/interval_pixels);
        imstack = zeros(nrow, ncolumn, ncoronal, 'uint16');
        for k = 1:nhorizontal
            try
                im = imread([channel_path filesep image_names{k}], 'tiff');
                imstack(k, :, :) = im(1:interval_pixels:end, :)';
            end
            disp([channel_path filesep image_names{k}]);
        end
        for k = 1:ncoronal
            im = imstack(:, :, k);
            save([save_path{i} filesep 'Channel' num2str(j) '_' num2str(k, '%03d') '.mat'], 'im');
        end
    end

    parfor j = 1:ncoronal
        for k = 1:nchannels
            im = load([save_path{i} filesep 'Channel' num2str(k) '_' num2str(j, '%03d') '.mat']);
            im = im.im;
            if k == 1
                imwrite(im, [save_path{i} filesep 'highres_' num2str(j, '%03d') '.tiff'], 'tiff', 'Compression', 'none');
            else
                imwrite(im, [save_path{i} filesep 'highres_' num2str(j, '%03d') '.tiff'], 'tiff', 'Compression', 'none', 'WriteMode', 'append');
            end
            delete([save_path{i} filesep 'Channel' num2str(k) '_' num2str(j, '%03d') '.mat']);
        end
    end

    workbar(i/numel(data_path), [num2str(i) '/' num2str(numel(data_path))], 'Progress');
end
msgbox('Done !');