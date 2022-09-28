upst_dir = '.\images'
downst_dir = '.\images\done'
images = dir(fullfile(upst_dir,'*.jpg'));

load('p4nonprocal') %from burnett et al 2018 supp info

for image = 1:length(images);
    if image == 1
        "now processing:"
    end
    images(image)
    [pathstr, name, ext] = fileparts(images(image).name)
    path = fullfile(upst_dir,strcat(name,ext));
    I = imread(path);
    [newI,conv]  = undistortImage(I,calibrationSession);
    imwrite(I,fullfile(downst_dir,strcat(name,ext)))
    imwrite(newI,fullfile(downst_dir,strcat(name,'_undist',ext)))
end
