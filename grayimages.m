dbstop if error

targetpath = './benchmark/targets/';
graypath = './benchmark/targets_gray/';

for i=1:500
    im = imread([targetpath, sprintf('%04d.jpg', i)]);
    imgray = rgb2gray(im);
    imwrite(imgray, [graypath, sprintf('gray_%04d.jpg', i)]);
end