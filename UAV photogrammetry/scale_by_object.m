function [range] = scale_by_object (image)
%crops image
%img1 = imrotate(img1,90,'bilinear');
imshow(image);
message = sprintf(['Bound the scaling object with a square and right click and click "crop image" when you are done.']);
uiwait(msgbox(message));
I2 = imcrop(this_frame_undist);
close all
img_wiener = imageenhancer(I2); %enhances image

%measurement
imshow(img_wiener);
message = sprintf(['Click a point on each end of the scaling object']);
test = msgbox(message);
mpos = get(test,'OuterPosition');
set(test,'OuterPosition',[mpos(1)+400,mpos(2:4)]);
uwait=(test);
hold on
n = 2
xlength=zeros(n,1);
ylength=zeros(n,1);
for ii = 1:n;
    [x1,y1]=ginput(1);
    plot(x1,y1,'--.r','MarkerSize',12);
    xlength(ii)=x1;
    ylength(ii)=y1;
end
tag_lenpx = sqrt((xlength(1)-xlength(2))^2+(ylength(1)-ylength(2))^2);
plot(xlength,ylength,'LineWidth',2);
dx=-5;dy=-5;
text(xlength+dx,ylength+dy,'OL');

% % plot measurements
% % save image code from: http://stackoverflow.com/questions/1848176/how-do-i-save-a-plotted-image-and-maintain-the-original-image-size-in-matlab
% [pathstr,name,ext] = fileparts(fullFileName);
% plottitle = [' Flight: ' num2str(flight) ' File: ' name ' Image #: ' num2str(wcnt)];
% title(plottitle)
% f = getframe(gcf);              %# Capture the current window
% imgout = [fullFileName(1:end-4),'_measured.jpg']
% imwrite(f.cdata,imgout);  %# Save the frame data

message = sprintf(['Press Enter twice to continue']);
test = msgbox(message);
pause
close all
    end