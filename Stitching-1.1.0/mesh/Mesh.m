classdef Mesh < handle
    %MESH Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        imgRows,imgCols;
        meshWidth,meshHeight;
        quadWidth,quadHeight;
        xMat,yMat;
    end
    
    methods 
        function obj = Mesh(rows,cols,quadWidth,quadHeight)
           obj.imgRows = rows;
           obj.imgCols = cols;
           obj.quadWidth = quadWidth;
           obj.quadHeight = quadHeight;
           
           xSet = 0;
           ySet = 0;
            
            x = 1;i = 1;
            while obj.imgCols - x > 0.5*quadWidth
                xSet(i) = x;
                i = i+1;
                x = x + obj.quadWidth;
            end
            xSet(i) = obj.imgCols;
            
            y =1; j = 1;
            while obj.imgRows - y >0.5*quadHeight
                ySet(j) = y;
                j = j + 1;
                y = y + obj.quadHeight;
            end
            ySet(j) = obj.imgRows;
            
            obj.meshWidth = length(xSet);
            obj.meshHeight =length(ySet);
            
            obj.xMat = zeros(obj.meshHeight,obj.meshWidth);
            obj.yMat = zeros(obj.meshHeight,obj.meshWidth);
            
            for y=1:obj.meshHeight
                for x=1:obj.meshWidth
                   obj.xMat(y,x) = xSet(x);
                   obj.yMat(y,x) = ySet(y);
                end
            end
        end
        function pt = getVertex(obj,i,j)
            pt = myPoint(obj.xMat(i+1,j+1),obj.yMat(i+1,j+1));
        end
        function setVertex(obj,i,j,pt)
           obj.xMat(i+1,j+1) = pt.x;
           obj.yMat(i+1,j+1) = pt.y;
        end
        function qd = getQuad(obj,i,j)
             V00 = obj.getVertex(i-1,j-1);
             V01 = obj.getVertex(i-1,j);
             V10 = obj.getVertex(i,j-1);
             V11 = obj.getVertex(i,j);
             qd = Quad(V00,V01,V10,V11);
        end
        function meshImg = drawMesh(obj,Img,gap,ab)
		  if ab == 'a'
		    color = [0 100 255];
		  else
		    color = [247 104 0];
		  end
          shapePoints = vision.ShapeInserter('Shape','Circles',...
                                     'Fill',true,...
                                     'FillColor','Custom',...
                                     'CustomFillColor',color,...
                                     'Opacity',0.8);
          shapeLine = vision.ShapeInserter('Shape','Lines',...
                                     'BorderColor','Custom',...
                                     'CustomBorderColor',color,...
                                     'Antialiasing',false);
          meshImg = Img;
          
          Points = zeros(obj.meshWidth*obj.meshHeight,3);
%           Lines = zeros(obj.meshWidth*obj.meshHeight,4);
          cc = 1;
          for i=1:obj.meshHeight
              for j=1:obj.meshWidth
                 Points(cc,1) = obj.xMat(i,j)+gap;
                 Points(cc,2) = obj.yMat(i,j)+gap;
                 Points(cc,3) = 5;
                 cc = cc + 1;
              end
          end
          
          cc=1;
          for i=1:obj.meshHeight-1
              for j=1:obj.meshWidth
                 Lines(cc,1) = obj.xMat(i,j)+gap;
                 Lines(cc,2) = obj.yMat(i,j)+gap;
                 Lines(cc,3) = obj.xMat(i+1,j)+gap;
                 Lines(cc,4) = obj.yMat(i+1,j)+gap;
                 cc = cc + 1;
              end
          end
          
          for i=1:obj.meshHeight
              for j=1:obj.meshWidth-1
                 Lines(cc,1) = obj.xMat(i,j)+gap;
                 Lines(cc,2) = obj.yMat(i,j)+gap;
                 Lines(cc,3) = obj.xMat(i,j+1)+gap;
                 Lines(cc,4) = obj.yMat(i,j+1)+gap;
                 cc = cc + 1;
              end
          end
          
          meshImg = step(shapeLine,meshImg,uint16(Lines));
          meshImg = step(shapePoints,meshImg,uint16(Points));
          
        end
    end
end

