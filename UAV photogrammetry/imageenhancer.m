function [img_wiener] = imageenhancer(I2)%enhances cropped image
I2 = im2double(I2);
srgb2lab = makecform('srgb2lab');
lab2srgb = makecform('lab2srgb');
img_lab = applycform(I2, srgb2lab); % convert to L*a*b*
% the values of luminosity can span a range from 0 to 100; scale them
% to [0 1] range (appropriate for MATLAB(R) intensity images of class double)
% before applying the three contrast enhancement techniques
max_luminosity = 100;
L = img_lab(:,:,1)/max_luminosity;
% replace the luminosity layer with the processed data and then convert
% the image back to the RGB colorspace
img_adapthisteq = img_lab;
img_adapthisteq(:,:,1) = adapthisteq(L)*max_luminosity;
img_adapthisteq = applycform(img_adapthisteq, lab2srgb);

%denoises image
img_wiener = img_adapthisteq;
img_wiener(:,:,1) = wiener2(img_adapthisteq(:,:,1),[3 3]);
img_wiener(:,:,2) = wiener2(img_adapthisteq(:,:,2),[3 3]);
img_wiener(:,:,3) = wiener2(img_adapthisteq(:,:,3),[3 3]);
end 