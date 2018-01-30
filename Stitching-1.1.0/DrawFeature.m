function I = DrawFeature( I , f, color, shape)
    I = insertMarker(I, f, shape, 'color', color, 'size', 2);        
    I = insertMarker(I, f, shape, 'color', color, 'size', 1);
    I = insertMarker(I, f, '*', 'color', color, 'size', 2);        
    I = insertMarker(I, f, '+', 'color', color, 'size', 2);        
    I = insertMarker(I, f, shape, 'color', color, 'size', 3);    
%     I = insertMarker(I, f, shape, 'color', color, 'size', 4);   
    I = insertMarker(I, f, shape, 'color', 'black', 'size', 4);
end

