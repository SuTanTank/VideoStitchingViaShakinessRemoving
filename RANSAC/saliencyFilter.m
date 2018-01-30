function [ f1, f2 ] = saliencyFilter(f1, f2, index )
%SALIENCYFILTER Summary of this function goes here
%   Detailed explanation goes here
    return
    filename = ['./SANY0025/' int2str(index) '_res.PNG'];
    sal = imread(filename);
    i = 1;
    if f2 ~= 0
        while i <= size(f1, 1)
            x1 = round(f1(i, 1)); 
            y1 = round(f1(i, 2)); 
            x2 = round(f2(i, 1)); 
            y2 = round(f2(i, 2)); 
            if (sal(y1, x1) > 70) || (sal(y2, x2) > 70)
                f1(i,:) = [];
                f2(i,:) = [];
            else
                i = i + 1;
            end
        end
    else
        while i <= size(f1, 1)
            x1 = round(f1(i, 1)); 
            y1 = round(f1(i, 2)); 
            if (sal(y1, x1) > 70)
                f1(i,:) = [];
            else
                i = i + 1;
            end
        end 
    end
end

