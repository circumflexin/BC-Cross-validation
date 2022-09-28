clear all

%% config
upstream_dir = 'D:\Raw\CAS 2019\Data\UAV\2019_09_08\1_0745'; % folder with your raw video, flight log and LiDAR data 
downstream_dir = 'D:\Analysis\Cross validation\UAV photogrammetry\outputs'; % folder you want your analysed images to go in
video_name = '05_52_a.MOV'; % make sure to check for additional chunks.
whale_id = "285"; % this is the 3s sighting number but any string (e.g. "gunslinger") will work. 
julian_day = "251";
UAV_flight_id = "UAV19_251_1"; % unique ID for cross referencing. Its also good to include this in you analysis filename e.g. UAV_AB_UAV19_251_1.m
tags = join(["SW19_250a",",","SW19_250b"]); % any associated tag records. Be consistent about whether you include these at the whale or deployment level.
UAV_logfile = 'DJIFlightRecord_2019-09-08_[07-53-19]-verb.csv'; %name of the flight log data (you'll need to convert these: https://www.phantomhelp.com/logviewer/upload/)
GPS_image = "JPG_05_52_b.JPG"; % image of the GPS enabled device for time sync
lens_cal = 'D:\Analysis\Cross validation\UAV photogrammetry\p4nonprocal'; % lens calibration object, you can calibrate the lens yourself in matlab (google image analysis toolbox) or use the supp. info from Torres et al 2018 
LiDAR_logfile = "UAV19_251_1_lid.csv"; % the CSV file from the lidar datalogger. You'll need to convert the time variables into a proper datetime format. See the attached example for an excel formula. 
LiDAR_timeout = duration(00,00,01);% only use LiDAR readings at least this close to the video frame. Can be increased in good conditions
%% required directories are created if they dont exist already
if ~exist(fullfile(downstream_dir,UAV_flight_id), 'dir')
    mkdir(fullfile(downstream_dir,UAV_flight_id))
end
if ~exist(fullfile(downstream_dir,UAV_flight_id,'/frames'), 'dir')
    mkdir(fullfile(downstream_dir,UAV_flight_id,'/frames'))
end
if ~exist(fullfile(downstream_dir,UAV_flight_id,'/corrected_frames'), 'dir')
    mkdir(fullfile(downstream_dir,UAV_flight_id,'/corrected_frames'))
end
if ~exist(fullfile(downstream_dir,UAV_flight_id,'/measurements'), 'dir')
    mkdir(fullfile(downstream_dir,UAV_flight_id,'/measurements'))
end
% check if the frame annotations file exists, create it if not
if ~isfile(fullfile(downstream_dir,UAV_flight_id,join(UAV_flight_id,'_frames.csv')))
    writetable(array2table(zeros(0,3),'VariableNames',{'Frame','Item','Note'}),fullfile(downstream_dir,UAV_flight_id,join([UAV_flight_id,"_frames.csv"])))
end

%% data collected by the UAV are loaded and cleaned
UAV = readtable(fullfile(upstream_dir,UAV_logfile));
IMU_time = table2array(UAV(:,1));
UAV_pitch = table2array(UAV(:,20));
UAV_roll = table2array(UAV(:,21));
UAV_tilt = rad2deg(atan(sqrt(tan(deg2rad(UAV_roll)).^2+tan(deg2rad(UAV_pitch)).^2))); %this is the one part of the analysis I am not 100% clear on. If you know anyone who could check the trigonometry I would be v grateful. 
tilt_bias = cos(deg2rad(UAV_tilt));


% only one datetime input format can be specificed at once, but round seconds are recorded without milliseconds, so are converted for consistency:
IMU_time_whole = datetime(IMU_time(:,1),'InputFormat','yyyy/MM/dd HH:mm:ss');
IMU_time_milli = datetime(IMU_time(:,1),'InputFormat','yyyy/MM/dd HH:mm:ss.SSS');
for i = 1 : length(IMU_time_milli)
    if isnat(IMU_time_milli(i,1))
        IMU_time_milli(i,1) = IMU_time_whole(i);
    end
end
IMU_time_milli = datetime(IMU_time_milli,'format','yyyy/MM/dd HH:mm:ss.SSS');
OSD_height = table2array(UAV(:,15));
plot(IMU_time_milli,OSD_height); % height:time as measured by the UAV. 

%  LiDAR data imported 
LiDAR = readtable(fullfile(upstream_dir,LiDAR_logfile),'DatetimeType','text');
lidar_time = datetime(LiDAR.Datetime(:,1),'InputFormat','dd/MM/yyyy HH:mm:ss');
lidar_height = LiDAR.Lidar_Height_Mtr
lidar_height(lidar_height >= 130) = NaN % readings of 130m are actually timeouts
plot(lidar_time,lidar_height) % inspect



% calculate and apply UAV:handheld GPS time offset:
im_info = imfinfo(fullfile(upstream_dir,GPS_image)) %exif data is extracted
imshow(fullfile(upstream_dir,GPS_image))
gps_time = datetime('08-Sep-2019 05:52:32'); % look at the image and input the time on the GPS here e.g 11-Sep-2019 07:43:12
offset = (datetime(im_info.DateTime, 'InputFormat', 'yyyy:MM:dd HH:mm:ss') - duration('02:00:00')) - gps_time % Offset is calculated. Make sure you account for timezone differences manually as here. 
IMU_time_milli_corr = IMU_time_milli - offset; % Drone time is synced with GPS time. 
for i = 1:length(lidar_time)
    [minValue, clostestindex] = min(abs(IMU_time_milli_corr-lidar_time(i)));
    if minValue < duration(00,00,01); 
        lidar_height_corr(i) = lidar_height(i) * tilt_bias(clostestindex); 
    else 
        lidar_height_corr(i) = NaN; 
    end
end

% and inspect, check the inflection point at takeoff or another rapid change in altitude to make sure the sync
% has worked
plot(lidar_time,lidar_height);
hold on
plot(IMU_time_milli_corr,OSD_height);
hold on
plot(lidar_time,lidar_height_corr);


% video file is loaded and metadata extracted
vid = VideoReader(fullfile(upstream_dir,video_name))

% Run through the video in slow motion, pause and note suitable frame numbers in the CSV in
% dataframe format. You should also find a frame with the hendheld GPS in it, ideally showing a second rollover.
% It's also good practice to score photogrammetry frames for suitabiliy
% following Christiansen et al 2016. If in doubt about a frame, include it. You
% can always drop it later. See the example for what codes to use for each
% type of frame (e.g. pgram for photogrammetry frame). 
implay(fullfile(upstream_dir,video_name),10) % video is played in slow motion
winopen(fullfile(downstream_dir,UAV_flight_id,join([UAV_flight_id,"_frames.csv"]))) % frame annotations file is opened in excel
vid_sync = table(datetime('08-Sep-2019 05:52:40'),33) % input the datetime and the frame number from the frame showing the GPS

% the CSV is loaded into MATLAB:
frame_tab = readtable(fullfile(downstream_dir,UAV_flight_id,join([UAV_flight_id,"_frames.csv"])));
frames = [] 
elaps_times = duration();
abs_times = datetime();
% alt_BAR = []
% range_LiD = []

load(lens_cal) 

i = 1;
for n = 1:length(frame_tab.Item) % for each frame in the CSV:
    if strcmp(frame_tab.Item(n), 'pgram') % this string signifies a good all round photogrammetry frame
        notes(i,1) = frame_tab.Note(n)
        frame_n = frame_tab.Frame(n);
        elaps_time = seconds(((frame_n)/vid.FrameRate)); % the elapsed time is calculated
        %elaps_time.Format = 'mm:ss.SSS';
        frames(i,1) = frame_n;
        elaps_times(i,1) = elaps_time;
        abs_time = datetime(vid_sync.Var1(1)) + (elaps_time -(seconds(vid_sync.Var2(1)/vid.FrameRate))); % the coordinated absolute time is calculated.
        abs_times(i,1) = abs_time;
        this_frame = read(vid,frame_n); % the frame frame is extracted 
        basename = join([UAV_flight_id, '_',frame_n,'.jpg'],""); % and saved
        name(i,1) = basename;
        imwrite(this_frame,fullfile(downstream_dir,UAV_flight_id,'/frames',basename),'jpg');
        [this_frame_undist,conv] = undistortImage(this_frame,calibrationSession); % distortion is corrected using coefficients from the lens calibration object. 
        imwrite(this_frame_undist,fullfile(downstream_dir,UAV_flight_id,'/corrected_frames',join([UAV_flight_id,'_',frame_n,'_','undist.jpg'],"")),'jpg');
        [minValue, clostestindex] = min(abs(IMU_time_milli_corr-abs_time)) % closest altitude is found
        frame_alt_BAR = OSD_height(clostestindex)
        alt_BAR(i,1) = frame_alt_BAR
        [minValue, clostestindex] = min(abs(lidar_time-abs_time))
        if minValue < LiDAR_timeout % if the frame is close enough in time to a lidar reading (specified in config section)
            range_LiD(i,1) = lidar_height(clostestindex) % the distance to the whale is calcualated
            range_LiD_corr(i,1) = lidar_height_corr(clostestindex)
        else 
            range_LiD(i,1) = NaN % otherwise, null values are returned. 
            range_LiD_corr(i,1) = NaN
        end
        abs_times.Format = 'yyyy/MM/dd HH:mm:ss.SSS'
        elaps_times.Format = 'mm:ss.SSS'
        % sum_tab = table([elaps_times,abs_times])
        % sense check using length of tag
        imshow(this_frame_undist)
        x = questdlg('Scaling object visble?') % if there's a scaling object in the frame you can use this to asess accuracy or calcuate the scaling factor independently. 
        close all
        switch x
            case 'Yes'
                tag_lenpx = scale_by_object_AB(this_frame_undist)
                %%
                %focal length
                cfl = 3.61; %mm  change this if not using a Phantom 4
                %pixel pitch (mm):
                % 6.3 / vid.width  change this if not using a Phantom 4
                pmm = 0.001640625;% change this if not using a Phantom 4
                tag_len_cm = 28% insert length of scaling object
                GSD = tag_len_cm / tag_lenpx %cm/px
                %GSD = (alt/cfl)*pmm %from http://www.asprs.org/a/publications/pers/2004journal/march/2004_mar_297-300.pdf
                % so:
                range = ((cfl*GSD)/pmm)/100 % rearranged and solved for range to whale. in m
                tag_scaled_ranges(i,1) = range;
            case 'No'
                tag_scaled_ranges(i,1) = NaN;
        end
        i = i + 1
    end
end

close all
close all hidden

% summary table for further analysis is written

summary = table(frames,tag_scaled_ranges, alt_BAR, range_LiD,range_LiD_corr,notes)
summary.pixpitch(:) = pmm
summary.foc_len(:) = cfl
summary.tags_on(:) = tags
summary.whale_id(:) = whale_id
summary.flight(:) = UAV_flight_id


writetable(summary, fullfile(downstream_dir,UAV_flight_id,strcat(UAV_flight_id,"_summary.csv"))) % comment this when you're done to prevent accidental overwriting.


