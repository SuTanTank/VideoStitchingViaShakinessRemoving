video = VideoReader('./case17-l.mp4');
if ~exist('left', 'dir')
    mkdir('left');
end
k = 0;
while hasFrame(video)
    k = k + 1;
    frame = readFrame(video);
    filename = ['./left/' sprintf('%03d',k) '.jpg'];
    imwrite(frame, filename);
end
video = VideoReader('./case17-r.mp4');
if ~exist('right', 'dir')
    mkdir('right');
end
k = 0;
while hasFrame(video)
    k = k + 1;
    frame = readFrame(video);
    filename = ['./right/' sprintf('%03d',k) '.jpg'];
    imwrite(frame, filename);
end