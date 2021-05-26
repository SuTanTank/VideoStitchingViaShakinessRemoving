% Read A1 C1A B1 C1B
% Get seam area
% sift A1 B1, filter outliers
% compute stitching error

addpath('RANSAC');
addpath('peter');
data = 'D:/MATLABs/Stitching/case19/res_only_2/seam/';
head = 1;
tail = 500;
nFrames = tail - head + 1;
ppf = 1000;
avgErrors = zeros(nFrames, 1);
for t = head:tail
    disp(['Render ' int2str(t) ]);
    A = imread([data '../' 'A' int2str(t) '.jpg']);
    A = imresize(A, 0.5);
    B = imread([data '../' 'B' int2str(t) '.jpg']);
    B = imresize(B, 0.5);
    [H, W, ~] = size(A);
    Acut = imread([data './' 'C' int2str(t) 'A.jpg']);
    Acut = imresize(Acut, 0.5);
    Bcut = imread([data './' 'C' int2str(t) 'B.jpg']);
    Bcut = imresize(Bcut, 0.5);
%     [f1, f2] = SIFT(A, B);
    [f1, f2] = SURF(A, B);

    A = im2double(A);
    B = im2double(B);
    AMask = (A > 0.05);
    AMask = AMask(:,:,1) | AMask(:,:,2) | AMask(:,:,3);
    AMask = double(AMask);
    BMask = (B > 0.05);
    BMask = BMask(:,:,1) | BMask(:,:,2) | BMask(:,:,3);
    BMask = double(BMask);


    ACutMask = (Acut > 0.1);
    ACutMask = ACutMask(:,:,1) | ACutMask(:,:,2) | ACutMask(:,:,3);
    ACutMask = double(ACutMask);
    BCutMask = (Bcut > 0.1);
    BCutMask = BCutMask(:,:,1) | BCutMask(:,:,2) | BCutMask(:,:,3);
    BCutMask = double(BCutMask);


    ACutMask = imgaussfilt(ACutMask .* BMask, 10);
    BCutMask = imgaussfilt(BCutMask .* AMask, 10);
    mask = ACutMask > 0 & BCutMask > 0;
    mask = mask & AMask & BMask;
%     figure(1);imshow(ACutMask);
%     figure(2);imshow(BCutMask);
%     imshow(mask);

    valid = mask(round(f1(:, 1) - 1) * H + round(f1(:, 2))) & mask(round(f2(:, 1) - 1) * H + round(f2(:, 2)));
    nf = sum(valid);
    f1m = f1(valid, :);
    f2m = f2(valid, :);

%     IA = insertMarker(A, f1m, 'o', 'color', 'red');
    IA = DrawFeature(A, f1m, 'red', 's');
    IB = DrawFeature(B, f1m, 'green', 's');
%     IB = insertMarker(B, f2m, 's', 'color', 'yellow');
    imshow(IA);
%     imshow(IB);

    d = f1m - f2m;

    error = (d(:, 1).^2 + d(:, 2).^2).^(0.5);
    error = sort(error);
    avgErrors(t) = mean(error(1:floor(nf * 10/ 10)));
    disp(avgErrors(t));
end

mean(avgErrors(~isnan(avgErrors)))