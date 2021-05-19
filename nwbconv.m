%%
% NWB file generator command line tool

function nwbconv(varargin)

    % set version number
    versionNumber = '0.1';

    % add script path
%{
    if ~isdeployed % checking MATLAB mode or stand-alone mode.
        [st,ind] = dbstack('-completenames');
        relpath = st(ind).file;
        [exedir,exename,ext] = fileparts(relpath);
        if exist([exedir '/util'],'dir')
            addpath([exedir '/util']);
            addpath([exedir '/lib']);
        end
    end
%}
    % get exe file full path
    global exePath;
    global exeName;
    [exePath, exeName, ext] = exeFilename();

    % init command line input
    handles.commandError = 0;
    handles.xlsFiles = {};
    handles.outpath = '.';
    handles.session = '1';

    % load command line input
    i = 1;
    while true
        if i > size(varargin, 2)
            break;
        end
        switch varargin{i}
            case {'--outpath'}
                handles.outpath = varargin{i+1};
                i = i + 1;
            case {'--session'}
                handles.session = varargin{i+1};
                i = i + 1;
            case {'-h','--help'}
                showUsage();
                return;
            case {'-v','--version'}
                disp([exeName ' version : ' num2str(versionNumber)]);
                return;
            otherwise
                if strcmp(varargin{i}(1), '-')
                    disp(['bad option : ' varargin{i}]);
                    i = size(varargin, 2);
                    handles.commandError = 1;
                else
                    handles.xlsFiles = [handles.xlsFiles varargin{i}];
                end
        end
        i = i + 1;
    end
    
    % check command input
    if handles.commandError
        showUsage();
        return;
    elseif isempty(handles.xlsFiles)
        disp('no input files. please specify nwb registration xlsx files.');
        showUsage();
        return;
    end

    % process input files
    processInputFiles(handles);
end

%%
% show usage function
function showUsage()
    global exePath;
    global exeName;
    disp(['usage: ' exeName ' [options] filename.xlsx ...']);
    disp('  --outpath path      output files path (default:".")');
    disp('  --session num       output session range (i.e. "1","1:8")');
    disp('  -v, --version       show version number');
    disp('  -h, --help          show command line help');
end

%%
% process input files (mail rutine)
%
function processInputFiles(handles)
    global exePath;
    global exeName;

    % init
    GEN_ROWS = 7;     % general row length
    SBJ_ROWS = 6;     % subject row length
    DEVICE_ROWS = 10; % device row length
    PHOTON_ROWS = 4;  % two photon series row length
    ROI_ROWS = 5;     % ROI row length
    ROIRESP_ROWS = 6; % ROI Response Series row length
    BE_ROWS = 7;      % behavior row length
    TRIAL_ROWS = 6;   % trials row length (variable)

    N = length(handles.xlsFiles);
    
    % process each file
    for i = 1:N
        % init data
        C1 = {};

        % load node status signals csv or mat file
        fname = handles.xlsFiles{i};
        if ~exist(fname,'file')
            disp(['file is not found. ignoring : ' fname]);
            continue;
        end
        [path,name,ext] = fileparts(fname);
        if strcmp(ext,'.xls') || strcmp(ext,'.xlsx')
            C1 = readcell(fname);
        else
            disp(['bad file format. ignoring : ' fname]);
            continue;
        end
        rowNum = size(C1,1);

        r=1;
        if ~strcmp(C1{r,1},'General Information')
            disp(['bad file context(' num2str(r) ',' num2str(1) '). ignoring : ' fname]);
            continue;
        end
        
        % parse session parameter
        st = 1;
        en = 1;
        range = strsplit(handles.session,':');
        if length(range) > 1
            st = str2num(range{1});
            en = str2num(range{2});
        elseif length(range{1}) > 0
            st = str2num(range{1});
            en = str2num(range{1});
        end
        % parse field
        for ii=st:en
            % replace params
            C = C1;
            for j=1:size(C,1)
                if isa(C{j,2},'char')
                    C{j,2} = strrep(C{j,2},'%session%',num2str(ii));
                end
            end
            
            % Set up the NWB file
            r=3; c=2;
            % option check
            if (length(C{r+3,c})==1&&ismissing(C{r+3,c}))||isempty(C{r+3,c}), C{r+3,c}=''; end
            if (length(C{r+4,c})==1&&ismissing(C{r+4,c}))||isempty(C{r+4,c}), C{r+4,c}=''; end
            if (length(C{r+5,c})==1&&ismissing(C{r+5,c}))||isempty(C{r+5,c}), C{r+5,c}=''; end
            if (length(C{r+6,c})==1&&ismissing(C{r+6,c}))||isempty(C{r+6,c}), C{r+6,c}=''; end
            nwb = NwbFile( ...
                'session_description', C{r,c},...
                'identifier', C{r+1,c}, ...
                'session_start_time', C{r+2,c}, ...
                'general_experimenter', C{r+3,c}, ... % optional
                'general_session_id', C{r+4,c}, ... % optional
                'general_institution', C{r+5,c}, ... % optional
                'general_related_publications', C{r+6,c}); % optional
            r = r + GEN_ROWS;

            % Subject information
            if ~strcmp(C{r,1},'Subject Information')
                disp(['bad file context(' num2str(r) ',' num2str(1) '). ignoring : ' fname]);
                continue;
            end
            subject = types.core.Subject( ...
                'subject_id', C{r+1,c}, ...
                'age', C{r+2,c}, ...
                'description', C{r+3,c}, ...
                'species', C{r+4,c}, ...
                'sex', C{r+5,c});
            nwb.general_subject = subject;
            r = r + SBJ_ROWS;

            % Optical Physiology
            if ~strcmp(C{r,1},'Optical Physiology')
                disp(['bad file context(' num2str(r) ',' num2str(1) '). ignoring : ' fname]);
                continue;
            end
            optical_channel = types.core.OpticalChannel( ...
                'description', C{r+5,c}, ...
                'emission_lambda', C{r+6,c});
            device = types.core.Device( ...
                'description', C{r+8,c}, ...
                'manufacturer', C{r+9,c});
            nwb.general_devices.set('Device', device);

            imaging_plane = types.core.ImagingPlane( ...
                'optical_channel', optical_channel, ...
                'description', C{r+1,c}, ...
                'device', types.untyped.SoftLink('/general/devices/Device'), ...
                'excitation_lambda', C{r+2,c}, ...
                'indicator', C{r+3,c}, ...
                'location', C{r+4,c});
            nwb.general_optophysiology.set('imaging_plane', imaging_plane);
            r = r + DEVICE_ROWS;

            isdata=1;
            if (length(C{r,c})==1&&ismissing(C{r,c}))||isempty(C{r,c}), isdata=0; end
            if isdata
                data = parseReadFile(C{r,c},3);
                image_series = types.core.TwoPhotonSeries( ...
                    'imaging_plane', types.untyped.SoftLink('/general/optophysiology/imaging_plane'), ...
                    'starting_time', C{r+2,c}, ...
                    'starting_time_rate', C{r+3,c}, ...
                    'data', data, ...
                    'data_unit', C{r+1,c});
                nwb.acquisition.set('TwoPhotonSeries', image_series);
            end
            r = r + PHOTON_ROWS;

            % Regions of interest (ROIs)
            while true
                if r > size(C,1) || ~strcmp(C{r,1},'Regions of interest')
                    break;
                end
                isdata=1; additional=0;
                if (length(C{r+4,c})==1&&ismissing(C{r+4,c}))||isempty(C{r+4,c}), isdata=0; end
                if isdata
                    % check colnames
                    colnames = parseColnames(C{r+2,c});
                    colD = {};
                    for j=2:length(colnames)
                        colD{j} = parseReadFile(C{r+3+j,c},1);
                    end

                    % read ROI masks
                    image_mask = parseReadFile(C{r+4,c},3);
                    n_rois = size(image_mask, 3);
                    ids = 0:(n_rois-1);
                    additional = length(colnames) - 1;
                    plane_segmentation = types.core.PlaneSegmentation( ...
                        'colnames', colnames, ...
                        'description', C{r+1,c}, ...
                        'id', types.hdmf_common.ElementIdentifiers('data', ids), ...
                        'imaging_plane', types.untyped.SoftLink('/general/optophysiology/imaging_plane'));
                    plane_segmentation.image_mask = types.hdmf_common.VectorData( ...
                        'data', image_mask, 'description', 'image masks');
                    varargin = {};
                    for j=2:length(colnames)
                        varargin{end+1} = colnames{j};
                        varargin{end+1} = types.hdmf_common.VectorData('data', colD{j}, 'description', [colnames{j} ' for ROI']);
                    end
                    if ~isempty(varargin)
                        [plane_segmentation.vectordata, ivarargin] = types.util.parseConstrained(plane_segmentation, 'vectordata', 'types.hdmf_common.VectorData', varargin{:});
                    end

                    img_seg = types.core.ImageSegmentation();
                    img_seg.planesegmentation.set('PlaneSegmentation', plane_segmentation);
                    ophys_module = types.core.ProcessingModule('description', 'contains optical physiology data');
                    ophys_module.nwbdatainterface.set('ImageSegmentation', img_seg);
                    nwb.processing.set('ophys', ophys_module);
                end
                r = r + ROI_ROWS + additional;
            end

            % ROI : Response Series
            while true
                if r > size(C,1) || ~strcmp(C{r,1},'ROI : Response Series')
                    break;
                end
                isdata=1;
                if (length(C{r+2,c})==1&&ismissing(C{r+2,c}))||isempty(C{r+2,c}), isdata=0; end
                if isdata
                    series_data = parseReadFile(C{r+2,c},1);
                    plane_seg_object_view = types.untyped.ObjectView( ...
                        '/processing/ophys/ImageSegmentation/PlaneSegmentation');
                    roi_table_region = types.hdmf_common.DynamicTableRegion( ...
                        'table', plane_seg_object_view, ...
                        'description', 'all_rois', ...
                        'data', [0 n_rois-1]);
                    roi_response_series = types.core.RoiResponseSeries( ...
                        'description', C{r+1,c}, ...
                        'rois', roi_table_region, ...
                        'data', series_data, ...
                        'data_unit', C{r+3,c}, ...
                        'starting_time', C{r+4,c}, ...
                        'starting_time_rate', C{r+5,c});
                    fluorescence = types.core.Fluorescence();
                    fluorescence.roiresponseseries.set('RoiResponseSeries', roi_response_series);
                    ophys_module.nwbdatainterface.set('Fluorescence', fluorescence);
                    nwb.processing.set('ophys', ophys_module);
                end
                r = r + ROIRESP_ROWS;
            end

            % Behavior : Spatial Series
            while true
                if r > size(C,1) || ~strcmp(C{r,1},'Behavior : Spatial Series')
                    break;
                end
                % check data cell
                isdata=1;
                if (length(C{r+2,c})==1&&ismissing(C{r+2,c}))||isempty(C{r+2,c}), isdata=0; end
                if isdata
                    position_data = parseReadFile(C{r+2,c},1);
                    spatial_series_ts = types.core.SpatialSeries( ...
                        'description', C{r+1,c}, ...
                        'data', position_data, ...
                        'data_unit', C{r+3,c}, ...
                        'starting_time', C{r+4,c}, ...
                        'starting_time_rate', C{r+5,c}, ...
                        'reference_frame', C{r+6,c});
                    Position = types.core.Position('SpatialSeries', spatial_series_ts);
                    % create processing module
                    behavior_mod = types.core.ProcessingModule( 'description',  'contains behavioral data');
                    % add the Position object (that holds the SpatialSeries object)
                    behavior_mod.nwbdatainterface.set('Position', Position);
                    % add the processing module to the NWBFile object, and name it "behavior"
                    nwb.processing.set('behavior', behavior_mod);
                end
                r = r + BE_ROWS;
            end

            % Trials
            while true
                if r > size(C,1) || ~strcmp(C{r,1},'Trials')
                    break;
                end
                % check id cell
                isdata=1; additional=0;
                if (length(C{r+4,c})==1&&ismissing(C{r+4,c}))||isempty(C{r+4,c}), isdata=0; end
                if isdata
                    % check colnames
                    colnames = parseColnames(C{r+2,c});
                    colD = {};
                    for j=3:length(colnames)
                        colD{j} = parseReadFile(C{r+3+j,c},1);
                    end
                    % check convert sequence to start&end
                    if contains(C{r+4,c},'#time_rate')
                        % get time rate
                        R = strsplit(C{r+4,c},'=');
                        rate = str2num(R{2});
                        % sequence data convert mode
                        [start_time, stop_time, colD] = convertSequenceToStartStop(colD, rate);
                    else
                        % read times
                        start_time = parseReadFile(C{r+4,c},1);
                        stop_time = parseReadFile(C{r+5,c},1);
                    end
                    ids = 0:(size(start_time)-1);
                    additional = length(colnames) - 2;
                    trials = types.core.TimeIntervals( ...
                        'colnames', colnames, ...
                        'description', C{r+1,c}, ...
                        'id', types.hdmf_common.ElementIdentifiers('data', ids), ...
                        'start_time', types.hdmf_common.VectorData('data', start_time, ...
                            'description','start time of trial'), ...
                        'stop_time', types.hdmf_common.VectorData('data', stop_time, ...
                            'description','end of each trial'));
                    varargin = {};
                    for j=3:length(colnames)
                        varargin{end+1} = colnames{j};
                        varargin{end+1} = types.hdmf_common.VectorData('data', colD{j}, 'description', [colnames{j} ' of each trial']);
                    end
                    if ~isempty(varargin)
                        [trials.vectordata, ivarargin] = types.util.parseConstrained(trials,'vectordata', 'types.hdmf_common.VectorData', varargin{:});
                    end
                    nwb.intervals_trials = trials;
                end
                r = r + TRIAL_ROWS + additional;
            end

            % export nwb file
            if st == 1 && en == 1
                session = '';
            else
                session = ['_' num2str(ii)];
            end
            fname = [handles.outpath '/' name session '.nwb'];
            disp(['export NWB file : ' fname]);
            nwbExport(nwb, fname);
        end
    end
end

%%
function [start_time, stop_time, colDout] = convertSequenceToStartStop(colD, rate)
    % init
    cueseries = colD{3};
    start_time = [];
    stop_time = [];
    colDout = {};
    for j=3:length(colD)
        colDout{j} = [];
    end

    % convert
    lastcue = 0; % cueseries(1); TODO:
    for i=1:length(cueseries)
        if cueseries(i) ~= lastcue
            if length(start_time) > length(stop_time)
                stop_time = [stop_time; (i-1) / rate];
                for j=3:length(colD)
                    A = colD{j}(start_i:i);
                    colDout{j} = [colDout{j}; mode(A)];
                end
            else
                start_time = [start_time; (i-1) / rate];
                start_i = i;
            end
            lastcue = cueseries(i);
        end
    end
    if length(start_time) > length(stop_time)
        stop_time = [stop_time; i / rate];
        for j=3:length(colD)
            A = colD{j}(start_i:i);
            colDout{j} = [colDout{j}; mode(A)];
        end
    end
    disp(['converting cue time series with time rate : ' num2str(rate)]);
end

%%
function C = parseColnames(str)
    C = strsplit(str,',');
end

%%
function data = parseReadFile(str, catDim)
    data = [];
    fileName = [];
    transpose = 0;
    logicalTh = [];
    remove = [];

    % command extraction
    command = strsplit(str,'#');
    for j=2:length(command)
        if contains(command{j}, 'concatdim') 
            C = strsplit(command{j},'=');
            catDim = str2num(C{2});
        elseif contains(command{j}, 'transpose') 
            transpose = 1;
        elseif contains(command{j}, 'remove') 
            C = strsplit(command{j},'=');
            D = strsplit(C{2},',');
            dim = str2num(D{1});
            E = strsplit(D{2},':');
            st = str2num(E{1});
            if length(E)>1, ed=str2num(E{2}); else, ed=st; end
            remove = [dim, st, ed];
        elseif contains(command{j}, 'logical') 
            C = strsplit(command{j},'=');
            logicalTh = str2num(C{2});
        end
    end

    % parse file string
    Cm = strsplit(command{1},';'); % file separation
    for j=1:length(Cm)
        % separate filename and fieldName
        C = strsplit(Cm{j},':');
        if length(C) == 1
            fileName = C{1};
            fieldName = 'data';
        elseif length(C) == 2
            fileName = C{1};
            D = strsplit(C{2},'.');
            fieldName = D{1};
            itemNames = D(2:end);
        end
        % check .mat or else
        if contains(fileName,'.mat')
            f = load(fileName);
            d = parseFieldNames(f, fieldName, itemNames, catDim);
            data = cat(catDim, data, d);
        elseif contains(fileName,'.tiff') ||  contains(fileName,'.tif')
            tstack  = Tiff(fileName);
            [I,J] = size(tstack.read());
            K = length(imfinfo(fileName));
            d = uint16(zeros(I,J,K));
            d(:,:,1)  = tstack.read();
            for n = 2:K
                tstack.nextDirectory()
                d(:,:,n) = tstack.read();
            end
            data = cat(catDim, data, d);
        elseif contains(fileName,'.png')
            d = imread(fileName);
            data = cat(catDim, data, d);
        elseif contains(fileName,'.csv')
            d = readmatrix(fileName);
            data = cat(catDim, data, d);
        end
    end
    % other operation
    if ~isempty(remove)
        switch remove(1)
        case 1
            if ndims(data)==2, data(remove(2):remove(3),:) = []; end
            if ndims(data)==3, data(remove(2):remove(3),:,:) = []; end
        case 2
            if ndims(data)==2, data(:,remove(2):remove(3)) = []; end
            if ndims(data)==3, data(:,remove(2):remove(3),:) = []; end
        case 3
            if ndims(data)==3, data(:,:,remove(2):remove(3)) = []; end
        end
    end
    if transpose
        order = [2 1];
        for j=3:length(size(data)), order=[order j]; end
        data = permute(data,order);
    end
    if ~isempty(logicalTh)
        data(data<logicalTh) = 0;
        data = logical(data);
    end

    % check data
    if isempty(data)
        disp(['failed to load data : ' str]);
    else
        disp(['success to load data : ' str]);
    end
end

%%
function data = parseFieldNames(f, fieldName, itemNames, catDim)
    data = [];
    % check cell or array
    if contains(fieldName, '(') || contains(fieldName, '{') 
        if contains(fieldName, '{'), iscell=1; else iscell=0; end
        if iscell
            F = strsplit(fieldName,'{');
            N = strsplit(F{2},'}');
        else
            F = strsplit(fieldName,'(');
            N = strsplit(F{2},')');
        end
        % parse number
        dim1 = []; dim2 = [];
        st = 1;
        en = length(f.(F{1}));
        R = strsplit(N{1},'-');
        if length(R) > 1
            st = str2num(R{1});
            en = str2num(R{2});
        elseif length(R{1}) > 0
            if contains(R{1}, ',')
                D = strsplit(R{1},',');
                dim1 = str2num(D{1});
                dim2 = str2num(D{2});
            else
                st = str2num(R{1});
                en = str2num(R{1});
            end
        end
        % parse field
        if isempty(dim1) && isempty(dim2)
            for ii=st:en
                if iscell
                    f2 = f.(F{1}){ii};
                else
                    f2 = f.(F{1})(ii);
                end
                if isempty(itemNames)
                    d = f2;
                else
                    d = parseFieldNames(f2, itemNames{1}, itemNames(2:end), catDim);
                end
                data = cat(catDim, data, d);
            end
        else
            % should be matrix data
            if isempty(dim1)
                d = f.(F{1})(:,dim2);
            elseif isempty(dim2)
                d = f.(F{1})(dim1,:);
            else
                d = f.(F{1})(dim1,dim2);
            end
            data = d;
        end
    else
        if isempty(itemNames)
            if contains(fieldName, ',')
                T = [];
                S = strsplit(fieldName,','); 
                for k=1:length(S)
                    tt = f.(S{k});
                    T = cat(2, T, tt(:));
                end
                d = T;
            else
                d = f.(fieldName);
                sz = size(d);
                if length(sz) > 2
                    d = squeeze(d);
                elseif sz(1) == 1 || sz(2) == 1
                    d = d(:);
                end
            end
        else
            f2 = f.(fieldName);
            d = parseFieldNames(f2, itemNames{1}, itemNames(2:end), catDim);
        end
        data = d;
    end
    % shrink data size
    if isa(data,'double')
        data = single(data);
    end
end

%%
function [exedir, exename, ext] = exeFilename()
    % Determine full path to calling file
    % In a deployed program, this will point to the application exe.
    % In regular Matlab, it will be the calling m-file.

    if isdeployed % Executable in stand-alone mode.
        [~, result] = system('path');
        exedir = char(regexpi(result, 'Path=(.*?);', 'tokens', 'once'));

        args = System.Environment.GetCommandLineArgs;
        relpath = char(args(1));
        [~,exename,ext] = fileparts(relpath);

    else % MATLAB mode.
        [st,ind] = dbstack('-completenames');
        relpath = st(ind+1).file;
        [exedir,exename,ext] = fileparts(relpath);
    end
end
