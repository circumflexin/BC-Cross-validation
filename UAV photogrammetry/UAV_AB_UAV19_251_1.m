
clear all

%% config
upstream_dir = 'D:\Raw\CAS 2019\Data\UAV\2019_09_08\1_0745';
downstream_dir = 'D:\Analysis\Cross validation\UAV photogrammetry\outputs';
video_name = '05_52_a.MOV';
whale_id = "285";
julian_day = "251";
UAV_flight_id = "UAV19_251_1";
tags = join(["SW19_250a",",","SW19_250b"]);
UAV_logfile = 'DJIFlightRecord_2019-09-08_[07-53-19]-verb.csv';
GPS_image = "JPG_05_52_b.JPG";
lens_cal = 'D:\Analysis\Cross validation\UAV photogrammetry\p4nonprocal';
LiDAR_logfile = "UAV19_251_1_lid.csv";
LiDAR_timeout = duration(00,00,01);% only use LiDAR readings at least this close to the video frame. Can be increased in good conditions
%% check if the required directories exist, create them if not

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

%% load and sanitise data collected onboard UAV
UAV = readtable(fullfile(upstream_dir,UAV_logfile));
IMU_time = table2array(UAV(:,1));
UAV_pitch = table2array(UAV(:,20));
UAV_roll = table2array(UAV(:,21));
UAV_tilt = rad2deg(atan(sqrt(tan(deg2rad(UAV_roll)).^2+tan(deg2rad(UAV_pitch)).^2)));
tilt_bias = cos(deg2rad(UAV_tilt));


% only one datetime input format can be specificed at once, but round seconds are recorded witout milliseconds, so:
IMU_time_whole = datetime(IMU_time(:,1),'InputFormat','yyyy/MM/dd HH:mm:ss');
IMU_time_milli = datetime(IMU_time(:,1),'InputFormat','yyyy/MM/dd HH:mm:ss.SSS');
for i = 1 : length(IMU_time_milli)
    if isnat(IMU_time_milli(i,1))
        IMU_time_milli(i,1) = IMU_time_whole(i);
    end
end
IMU_time_milli = datetime(IMU_time_milli,'format','yyyy/MM/dd HH:mm:ss.SSS');
OSD_height = table2array(UAV(:,15));
plot(IMU_time_milli,OSD_height);

% import and sanitise LiDAR data
LiDAR = readtable(fullfile(upstream_dir,LiDAR_logfile),'DatetimeType','text');
lidar_time = datetime(LiDAR.Datetime(:,1),'InputFormat','dd/MM/yyyy HH:mm:ss');
lidar_height = LiDAR.Lidar_Height_Mtr
lidar_height(lidar_height >= 130) = NaN % readings of 130m are actually timeouts
plot(lidar_time,lidar_height) % inspect



% calculate and apply UAV:handheld GPS time offset:
im_info = imfinfo(fullfile(upstream_dir,GPS_image))
imshow(fullfile(upstream_dir,GPS_image))
gps_time = datetime('08-Sep-2019 05:52:32'); % look at the image and input the time on the hh GPS e.g 11-Sep-2019 07:43:12
offset = (datetime(im_info.DateTime, 'InputFormat', 'yyyy:MM:dd HH:mm:ss') - duration('02:00:00')) - gps_time
IMU_time_milli_corr = IMU_time_milli - offset; % apply first round sync, accurate to nearest second
for i = 1:length(lidar_time)
    [minValue, clostestindex] = min(abs(IMU_time_milli_corr-lidar_time(i)));
    if minValue < duration(00,00,01);
        lidar_height_corr(i) = lidar_height(i) * tilt_bias(clostestindex);
    else 
        lidar_height_corr(i) = NaN;
    end
end

% and inspect, check the inflection point at takeoff to make sure the sync
% has worked
plot(lidar_time,lidar_height);
hold on
plot(IMU_time_milli_corr,OSD_height);
hold on
plot(lidar_time,lidar_height_corr);

% cross correlation
OSD_height_downs = downsample(OSD_height, 5);
X1 = xcorr(OSD_height_downs,lidar_height);

% [m,d]=max(X1);      %find value and index of maximum value of cross-correlation amplitude
% delay=d-max(length(OSD_height_downs),length(lidar_height));   %shift index d, as length(X1)=2*N-1; where N is the length of the signals
% figure,plot(length(OSD_height_downs))                                    %Plot signal s1
% hold on,plot([delay+1:length(OSD_height_downs)+delay],lidar_height,'r');   %Delay signal s2 by delay in order to align them
% grid on

% load video file
vid = VideoReader(fullfile(upstream_dir,video_name))

% run through the video in slow motion, pause and note suitable frame numbers in .CSV in
% dataframe format.
implay(fullfile(upstream_dir,video_name),10) % play the video in slow motion
winopen(fullfile(downstream_dir,UAV_flight_id,join([UAV_flight_id,"_frames.csv"]))) % tell excel to open the frame annotations file
vid_sync = table(datetime('08-Sep-2019 05:52:40'),33) %datetime in frame, frame number

%load the CSV
frame_tab = readtable(fullfile(downstream_dir,UAV_flight_id,join([UAV_flight_id,"_frames.csv"])));
frames = [] % to do: preallocate
elaps_times = duration();
abs_times = datetime();
% name = [] % to do: preallocate
% alt_BAR = []
% range_LiD = []

load(lens_cal) % load lens calibration file

i = 1;
for n = 1:length(frame_tab.Item)
    if strcmp(frame_tab.Item(n), 'pgram') % this string signifies a good all round photogrammetry frame
        notes(i,1) = frame_tab.Note(n)
        frame_n = frame_tab.Frame(n);
        elaps_time = seconds(((frame_n)/vid.FrameRate)); % from start of recording. to do: softcode this.
        %elaps_time.Format = 'mm:ss.SSS';
        frames(i,1) = frame_n;
        elaps_times(i,1) = elaps_time;
        abs_time = datetime(vid_sync.Var1(1)) + (elaps_time -(seconds(vid_sync.Var2(1)/vid.FrameRate)));
        abs_times(i,1) = abs_time;
        this_frame = read(vid,frame_n); % extract the frame
        basename = join([UAV_flight_id, '_',frame_n,'.jpg'],"");
        name(i,1) = basename;
        imwrite(this_frame,fullfile(downstream_dir,UAV_flight_id,'/frames',basename),'jpg');
        [this_frame_undist,conv] = undistortImage(this_frame,calibrationSession); % undistort using calibration coefficients from Burnett et al 2018
        imwrite(this_frame_undist,fullfile(downstream_dir,UAV_flight_id,'/corrected_frames',join([UAV_flight_id,'_',frame_n,'_','undist.jpg'],"")),'jpg');
        [minValue, clostestindex] = min(abs(IMU_time_milli_corr-abs_time)) % find closest altitude
        frame_alt_BAR = OSD_height(clostestindex)
        alt_BAR(i,1) = frame_alt_BAR
        [minValue, clostestindex] = min(abs(lidar_time-abs_time))
        if minValue < LiDAR_timeout
            range_LiD(i,1) = lidar_height(clostestindex)
            range_LiD_corr(i,1) = lidar_height_corr(clostestindex)
        else
            range_LiD(i,1) = NaN
            range_LiD_corr(i,1) = NaN
        end
        abs_times.Format = 'yyyy/MM/dd HH:mm:ss.SSS'
        elaps_times.Format = 'mm:ss.SSS'
        % sum_tab = table([elaps_times,abs_times])
        % sense check using length of tag
        imshow(this_frame_undist)
        x = questdlg('Tag visble?')
        close all
        switch x
            case 'Yes'
                tag_lenpx = scale_by_object_AB(this_frame_undist)
                %%
                %focal length
                cfl = 3.61; %mm
                %pixel pitch %mm
                % 6.3 / vid.width should this be diagonal?
                pmm = 0.001640625;% change to calculate from the image?
                tag_len_cm = 28% 28 with aerial, 24.5 without. 
                GSD = tag_len_cm / tag_lenpx %cm/px
                %GSD = (alt/cfl)*pmm %from http://www.asprs.org/a/publications/pers/2004journal/march/2004_mar_297-300.pdf
                % so:
                range = ((cfl*GSD)/pmm)/100 % rearrange and solve for range to whale. in m
                tag_scaled_ranges(i,1) = range;
            case 'No'
                tag_scaled_ranges(i,1) = NaN;
        end
        i = i + 1
    end
end

close all
close all hidden

summary = table(frames,tag_scaled_ranges, alt_BAR, range_LiD,range_LiD_corr,notes)
summary.pixpitch(:) = pmm
summary.foc_len(:) = cfl
summary.tags_on(:) = tags
summary.whale_id(:) = whale_id
summary.flight(:) = UAV_flight_id


%writetable(summary, fullfile(downstream_dir,UAV_flight_id,strcat(UAV_flight_id,"_summary.csv"))) % comment this when you're done to prevent overwriting.

