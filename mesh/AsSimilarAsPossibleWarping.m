classdef AsSimilarAsPossibleWarping < handle
    
    properties
        %smooth constraints
        SmoothConstraints;
        SCc;
        num_smooth_cons;
        
        %data constraints
        dataterm_element_i;
        dataterm_element_j;
        dataterm_element_orgPt;
        dataterm_element_desPt;
        dataterm_element_V00;
        dataterm_element_V01;
        dataterm_element_V10;
        dataterm_element_V11;
        DataConstraints;
        DCc;
        num_data_cons;
        
        
        rowCount;
        columns;
        
        x_index;
        y_index;
        
        %mesh
        source
        destin
        
        %control points
        sourcePts;
        targetPts;
        
        %histogram of the control points in each mesh
        histogram;
        threshold;
        weight;
        ADAPTIVE_WEIGHT;
        
        height,width; %mesh height,mesh width
        quadWidth,quadHeight; %quadWidth,%quadHeight
        imgHeight,imgWidth; %imgHeight,imgWidth
        
        warpIm;
        gap;
    end
    
    
    methods (Access = public)
        function obj = AsSimilarAsPossibleWarping(height,width,quadWidth,quadHeight,weight)
            obj.source = Mesh(height,width,quadWidth,quadHeight);
            obj.destin = Mesh(height,width,quadWidth,quadHeight);
            
            obj.imgHeight = height;
            obj.imgWidth = width;
            obj.quadWidth = quadWidth;
            obj.quadHeight = quadHeight;
            obj.height = obj.source.meshHeight;
            obj.width = obj.source.meshWidth;
            
            obj.threshold = quadWidth * quadHeight / 300;
            obj.weight = weight;
            
            obj.x_index = 0:obj.height*obj.width-1;
            obj.y_index = obj.height*obj.width:2*obj.height*obj.width-1;
            
            obj.num_smooth_cons = (obj.height-2)*(obj.width-2)*16 + (2*(obj.width+obj.height)-8)*8+4*4;% Total number of trangles * 2 (x and y)
            obj.columns = obj.width*obj.height*2;
            
            obj.SmoothConstraints = zeros(obj.num_smooth_cons,3);
            obj.SCc = 1;
            
            % obj.CreateSmoothCons(weight);
            % Move to setcontrolpts()
            
        end
        function SetControlPts(obj,inputsPts,outputsPts)
            len = length(inputsPts);
            
            obj.dataterm_element_orgPt = inputsPts;
            obj.dataterm_element_desPt = outputsPts;
            
            obj.dataterm_element_i = zeros(len,1);
            obj.dataterm_element_j = zeros(len,1);
            
            obj.dataterm_element_V00 = zeros(len,1);
            obj.dataterm_element_V01 = zeros(len,1);
            obj.dataterm_element_V10 = zeros(len,1);
            obj.dataterm_element_V11 = zeros(len,1);
            
            obj.histogram = zeros(obj.height - 1, obj.width - 1);
            
            for i=1:len
                pt = myPoint(inputsPts(i,1),inputsPts(i,2));
                obj.dataterm_element_i(i) = floor(pt.y/obj.quadHeight)+1;
                obj.dataterm_element_j(i) = floor(pt.x/obj.quadWidth)+1;
                % Find out which quad the point belongs to
                
                qd = obj.source.getQuad(obj.dataterm_element_i(i),obj.dataterm_element_j(i));
                % locate the quad
                
                coefficients = qd.getBilinearCoordinates(pt);
                obj.dataterm_element_V00(i) = coefficients(1);
                obj.dataterm_element_V01(i) = coefficients(2);
                obj.dataterm_element_V10(i) = coefficients(3);
                obj.dataterm_element_V11(i) = coefficients(4);
                % Calculate the bilinear coordinate
                
                obj.histogram(obj.dataterm_element_i(i),obj.dataterm_element_j(i)) = obj.histogram(obj.dataterm_element_i(i),obj.dataterm_element_j(i)) + 1;
            end
        end
        function Solve(obj)
            obj.CreateSmoothCons(obj.weight);
            b = CreateDataCons(obj);
            %              bb = b(obj.num_smooth_cons+1:end);
            N = length(obj.SmoothConstraints) + length(obj.DataConstraints);
            
            ARows = zeros(N,1);
            ACols = zeros(N,1);
            AVals = zeros(N,1);
            %              AARows = zeros(length(obj.DataConstraints), 1);
            %              AACols = zeros(length(obj.DataConstraints), 1);
            %              AAVals = zeros(length(obj.DataConstraints), 1);
            
            cc=1;
            for i=1:length(obj.SmoothConstraints)
                ARows(cc) = obj.SmoothConstraints(i,1)+1;
                ACols(cc) = obj.SmoothConstraints(i,2)+1;
                AVals(cc) = obj.SmoothConstraints(i,3);
                cc = cc + 1;
            end
            
            for i=1:length(obj.DataConstraints)
                ARows(cc) = obj.DataConstraints(i,1)+1;
                ACols(cc) = obj.DataConstraints(i,2)+1;
                AVals(cc) = obj.DataConstraints(i,3);
                cc = cc + 1;
            end
            
            A = sparse(ARows,ACols,AVals);
            
            x = A\b;
            
            
            halfcolumns = obj.columns/2;
            for i=0:obj.height-1
                for j=0:obj.width-1
                    pt = myPoint(x(i*obj.width+j+1),x(halfcolumns+i*obj.width+j+1));
                    obj.destin.setVertex(i,j,pt);
                end
            end
        end
        function warp = Warp(obj,Img,gap)
            
            obj.gap = gap;
            
            Img = double(Img);
            obj.warpIm = zeros(obj.imgHeight+gap*2,obj.imgWidth+gap*2,3);
            
            for i=1:obj.height-1
                for j=1:obj.width-1
                    p0 = obj.source.getVertex(i-1,j-1);
                    p1 = obj.source.getVertex(i-1,j);
                    p2 = obj.source.getVertex(i,j-1);
                    p3 = obj.source.getVertex(i,j);
                    
                    q0 = obj.destin.getVertex(i-1,j-1);
                    q1 = obj.destin.getVertex(i-1,j);
                    q2 = obj.destin.getVertex(i,j-1);
                    q3 = obj.destin.getVertex(i,j);
                    
                    qd1 = Quad(p0,p1,p2,p3);
                    qd2 = Quad(q0,q1,q2,q3);
                    quadWarp(obj,Img,qd1,qd2);
                end
            end
            
            warp = obj.warpIm;
        end
        function homos = CalcHomos(obj)
            homos = zeros(obj.height-1,obj.width-1,3,3);
            
            for i=1:obj.height-1
                for j=1:obj.width-1
                    q1 = obj.source.getQuad(i,j);
                    q2 = obj.destin.getQuad(i,j);
                    
                    source = zeros(4,2);
                    target = zeros(4,2);
                    
                    source(1,1) = q1.V00.x;source(1,2) = q1.V00.y;
                    source(2,1) = q1.V01.x;source(2,2) = q1.V01.y;
                    source(3,1) = q1.V10.x;source(3,2) = q1.V10.y;
                    source(4,1) = q1.V11.x;source(4,2) = q1.V11.y;
                    
                    target(1,1) = q2.V00.x;target(1,2) = q2.V00.y;
                    target(2,1) = q2.V01.x;target(2,2) = q2.V01.y;
                    target(3,1) = q2.V10.x;target(3,2) = q2.V10.y;
                    target(4,1) = q2.V11.x;target(4,2) = q2.V11.y;
                    
                    H = homography_4pts(source',target');
                    H = H./H(3,3);
                    
                    homos(i,j,:,:) = H;
                end
            end
            
        end
        function error = CalcError(obj)
            % e_h
            lens = length(obj.dataterm_element_desPt);
            error_h = 0;
            for i = 1:lens
                qd = obj.destin.getQuad(obj.dataterm_element_i(i),obj.dataterm_element_j(i));
                term_x = qd.V00.x * obj.dataterm_element_V00(i) ...
                    + qd.V01.x * obj.dataterm_element_V01(i) ...
                    + qd.V10.x * obj.dataterm_element_V10(i) ...
                    + qd.V11.x * obj.dataterm_element_V11(i) ...
                    - obj.dataterm_element_desPt(i,1);
                term_y = qd.V00.y * obj.dataterm_element_V00(i) ...
                    + qd.V01.y * obj.dataterm_element_V01(i) ...
                    + qd.V10.y * obj.dataterm_element_V10(i) ...
                    + qd.V11.y * obj.dataterm_element_V11(i) ...
                    - obj.dataterm_element_desPt(i,2);
                error_h = error_h + sqrt(term_x*term_x + term_y*term_y);
            end
            error_h = error_h / lens;
            
            %e_s
            error = error_h;
            
        end
    end
    
    methods (Access = private) %functions of solving as similar as possible warping
        function [u,v] = getSmoothWeight(obj,V1,V2,V3)
            d1 = sqrt((V1.x - V2.x)*(V1.x - V2.x) + (V1.y - V2.y)*(V1.y - V2.y));
            d3 = sqrt((V2.x - V3.x)*(V2.x - V3.x) + (V2.y - V3.y)*(V2.y - V3.y));
            
            v21 = myPoint(V1.x-V2.x,V1.y-V2.y);
            v23 = myPoint(V3.x-V2.x,V3.y-V2.y);
            
            cosin = v21.x*v23.x + v21.y*v23.y;
            cosin = cosin/(d1*d3);
            
            u_dis = cosin*d1;
            u = u_dis/d3;
            
            v_dis = sqrt(d1*d1 - u_dis*u_dis);
            v = v_dis/d3;
        end
        function CreateSmoothCons(obj,weight)
            obj.rowCount = 0;
            i=0;j=0;
            alpha = obj.getAlpha(i+1,j+1, weight);
            obj.addCoefficient_5(i,j,alpha);
            obj.addCoefficient_6(i,j,alpha);
            
            i=0;j=obj.width-1;
            alpha = obj.getAlpha(i+1,j, weight);
            obj.addCoefficient_7(i,j,alpha);
            obj.addCoefficient_8(i,j,alpha);
            
            i=obj.height-1;j=0;
            alpha = obj.getAlpha(i,j+1, weight);
            obj.addCoefficient_3(i,j,alpha);
            obj.addCoefficient_4(i,j,alpha);
            
            i=obj.height-1;j=obj.width-1;
            alpha = obj.getAlpha(i,j, weight);
            obj.addCoefficient_1(i,j,alpha);
            obj.addCoefficient_2(i,j,alpha);
            
            i=0;
            for j=1:obj.width-2
                alpha = mean([obj.getAlpha(i+1,j+1, weight) obj.getAlpha(i+1,j, weight)]);
                obj.addCoefficient_5(i,j,alpha);
                obj.addCoefficient_6(i,j,alpha);
                obj.addCoefficient_7(i,j,alpha);
                obj.addCoefficient_8(i,j,alpha);
            end
            
            i=obj.height-1;
            for j=1:obj.width-2
                alpha = mean([obj.getAlpha(i,j+1, weight) obj.getAlpha(i,j, weight)]);
                obj.addCoefficient_1(i,j,alpha);
                obj.addCoefficient_2(i,j,alpha);
                obj.addCoefficient_3(i,j,alpha);
                obj.addCoefficient_4(i,j,alpha);
            end
            
            j=0;
            for i=1:obj.height-2
                alpha = mean([obj.getAlpha(i+1,j+1, weight) obj.getAlpha(i,j+1, weight)]);
                obj.addCoefficient_3(i,j,alpha);
                obj.addCoefficient_4(i,j,alpha);
                obj.addCoefficient_5(i,j,alpha);
                obj.addCoefficient_6(i,j,alpha);
            end
            
            j=obj.width-1;
            for i=1:obj.height-2
                alpha = mean([obj.getAlpha(i,j, weight) obj.getAlpha(i+1,j, weight)]);
                obj.addCoefficient_1(i,j,alpha);
                obj.addCoefficient_2(i,j,alpha);
                obj.addCoefficient_7(i,j,alpha);
                obj.addCoefficient_8(i,j,alpha);
            end
            
            for i=1:obj.height-2
                for j=1:obj.width-2
                    alpha = mean([obj.getAlpha(i+1,j+1, weight) obj.getAlpha(i+1,j, weight) obj.getAlpha(i,j+1, weight) obj.getAlpha(i,j, weight)]);
                    obj.addCoefficient_1(i,j,alpha);
                    obj.addCoefficient_2(i,j,alpha);
                    obj.addCoefficient_3(i,j,alpha);
                    obj.addCoefficient_4(i,j,alpha);
                    obj.addCoefficient_5(i,j,alpha);
                    obj.addCoefficient_6(i,j,alpha);
                    obj.addCoefficient_7(i,j,alpha);
                    obj.addCoefficient_8(i,j,alpha);
                end
            end
            
            
        end
        function b = CreateDataCons(obj)
            len = length(obj.dataterm_element_i);% number of track points
            obj.num_data_cons = len*2;
            b = zeros(obj.num_data_cons + obj.num_smooth_cons,1);
            
            obj.DCc = 1;
            
            for k=1:len
                i = obj.dataterm_element_i(k);
                j = obj.dataterm_element_j(k);
                % locate the quad
                
                v00 = obj.dataterm_element_V00(k);
                v01 = obj.dataterm_element_V01(k);
                v10 = obj.dataterm_element_V10(k);
                v11 = obj.dataterm_element_V11(k);
                % get the coordinate of the quad
                
                obj.DataConstraints(obj.DCc,1) = obj.rowCount; obj.DataConstraints(obj.DCc,2) = obj.x_index((i-1)*obj.width+j-1+1);obj.DataConstraints(obj.DCc,3)= v00; obj.DCc = obj.DCc+1;
                obj.DataConstraints(obj.DCc,1) = obj.rowCount; obj.DataConstraints(obj.DCc,2) = obj.x_index((i-1)*obj.width+j+1);obj.DataConstraints(obj.DCc,3)= v01;   obj.DCc = obj.DCc+1;
                obj.DataConstraints(obj.DCc,1) = obj.rowCount; obj.DataConstraints(obj.DCc,2) = obj.x_index(i*obj.width+j-1+1);obj.DataConstraints(obj.DCc,3)= v10;     obj.DCc = obj.DCc+1;
                obj.DataConstraints(obj.DCc,1) = obj.rowCount; obj.DataConstraints(obj.DCc,2) = obj.x_index(i*obj.width+j+1);obj.DataConstraints(obj.DCc,3)= v11;       obj.DCc = obj.DCc+1;
                obj.rowCount = obj.rowCount+1;
                b(obj.rowCount) = obj.dataterm_element_desPt(k,1);
                
                obj.DataConstraints(obj.DCc,1) = obj.rowCount; obj.DataConstraints(obj.DCc,2) = obj.y_index((i-1)*obj.width+j-1+1);obj.DataConstraints(obj.DCc,3)= v00; obj.DCc = obj.DCc+1;
                obj.DataConstraints(obj.DCc,1) = obj.rowCount; obj.DataConstraints(obj.DCc,2) = obj.y_index((i-1)*obj.width+j+1);obj.DataConstraints(obj.DCc,3)= v01;   obj.DCc = obj.DCc+1;
                obj.DataConstraints(obj.DCc,1) = obj.rowCount; obj.DataConstraints(obj.DCc,2) = obj.y_index(i*obj.width+j-1+1);obj.DataConstraints(obj.DCc,3)= v10;     obj.DCc = obj.DCc+1;
                obj.DataConstraints(obj.DCc,1) = obj.rowCount; obj.DataConstraints(obj.DCc,2) = obj.y_index(i*obj.width+j+1);obj.DataConstraints(obj.DCc,3)= v11;       obj.DCc = obj.DCc+1;
                obj.rowCount = obj.rowCount+1;
                b(obj.rowCount) = obj.dataterm_element_desPt(k,2);
                
            end
        end
        function addCoefficient_1(obj,i,j,weight)
            %V3(i-1,j-1)
            %|
            %|____
            %V2   V1(i,j)
            V1 = obj.source.getVertex(i,j);
            V2 = obj.source.getVertex(i,j-1);
            V3 = obj.source.getVertex(i-1,j-1);
            
            [u,v] = obj.getSmoothWeight(V1,V2,V3);
            
            coordv1 =  i*obj.width+j;
            coordv2 =  i*obj.width+j-1;
            coordv3 = (i-1)*obj.width+j-1;
            
            coordv1 = coordv1 + 1;
            coordv2 = coordv2 + 1;
            coordv3 = coordv3 + 1;
            
            
            %            obj.SmoothCons.add_row(rowCount,obj.x_index(coordv2),(1.0-u)*weight);
            %            obj.SmoothCons.add_row(rowCount,obj.x_index(coordv3),u*weight);
            %            obj.SmoothCons.add_row(rowCount,obj.y_index(coordv2),v*weight);
            %            obj.SmoothCons.add_row(rowCount,obj.y_index(coordv3),(-1.0*v)*weight);
            %            obj.SmoothCons.add_row(rowCount,obj.x_index(coordv1),(-1.0)*weight);
            
            %V1.x =  V2.x + u * (V3.x - V2.x) + v * (V2.y - V3.y);
            obj.SmoothConstraints(obj.SCc,1) = obj.rowCount; obj.SmoothConstraints(obj.SCc,2) = obj.x_index(coordv2); obj.SmoothConstraints(obj.SCc,3) = (1.0-u)*weight;  obj.SCc = obj.SCc + 1;
            obj.SmoothConstraints(obj.SCc,1) = obj.rowCount; obj.SmoothConstraints(obj.SCc,2) = obj.x_index(coordv3); obj.SmoothConstraints(obj.SCc,3) = u*weight;        obj.SCc = obj.SCc + 1;
            obj.SmoothConstraints(obj.SCc,1) = obj.rowCount; obj.SmoothConstraints(obj.SCc,2) = obj.y_index(coordv2); obj.SmoothConstraints(obj.SCc,3) = v*weight;        obj.SCc = obj.SCc + 1;
            obj.SmoothConstraints(obj.SCc,1) = obj.rowCount; obj.SmoothConstraints(obj.SCc,2) = obj.y_index(coordv3); obj.SmoothConstraints(obj.SCc,3) = (-1.0*v)*weight; obj.SCc = obj.SCc + 1;
            obj.SmoothConstraints(obj.SCc,1) = obj.rowCount; obj.SmoothConstraints(obj.SCc,2) = obj.x_index(coordv1); obj.SmoothConstraints(obj.SCc,3) = (-1.0)*weight;   obj.SCc = obj.SCc + 1;
            obj.rowCount = obj.rowCount+1;
            
            
            %            obj.SmoothCons.add_row(rowCount,obj.y_index(coordv2),(1.0-u)*weight);
            %            obj.SmoothCons.add_row(rowCount,obj.y_index(coordv3),u*weight);
            %            obj.SmoothCons.add_row(rowCount,obj.x_index(coordv3),v*weight);
            %            obj.SmoothCons.add_row(rowCount,obj.x_index(coordv2),(-1.0*v)*weight);
            %            obj.SmoothCons.add_row(rowCount,obj.y_index(coordv1),(-1.0)*weight);
            %            rowCount = rowCount+1;
            
            %V1.y = V2.y + u * (V3.y - V2.y) + v * (V3.x - V2.x);
            obj.SmoothConstraints(obj.SCc,1) = obj.rowCount; obj.SmoothConstraints(obj.SCc,2) = obj.y_index(coordv2); obj.SmoothConstraints(obj.SCc,3) = (1.0-u)*weight;    obj.SCc = obj.SCc + 1;
            obj.SmoothConstraints(obj.SCc,1) = obj.rowCount; obj.SmoothConstraints(obj.SCc,2) = obj.y_index(coordv3); obj.SmoothConstraints(obj.SCc,3) = u*weight;          obj.SCc = obj.SCc + 1;
            obj.SmoothConstraints(obj.SCc,1) = obj.rowCount; obj.SmoothConstraints(obj.SCc,2) = obj.x_index(coordv3); obj.SmoothConstraints(obj.SCc,3) = v*weight;          obj.SCc = obj.SCc + 1;
            obj.SmoothConstraints(obj.SCc,1) = obj.rowCount; obj.SmoothConstraints(obj.SCc,2) = obj.x_index(coordv2); obj.SmoothConstraints(obj.SCc,3) = (-1.0*v)*weight;   obj.SCc = obj.SCc + 1;
            obj.SmoothConstraints(obj.SCc,1) = obj.rowCount; obj.SmoothConstraints(obj.SCc,2) = obj.y_index(coordv1); obj.SmoothConstraints(obj.SCc,3) = (-1.0)*weight;     obj.SCc = obj.SCc + 1;
            obj.rowCount = obj.rowCount+1;
            
            
            
        end
        function addCoefficient_2(obj,i,j,weight)
            %              V3   V2
            %             _____
            %                 |
            %                 |
            %                 V1(i,j)
            
            V1 = obj.source.getVertex(i,j);
            V2 = obj.source.getVertex(i-1,j);
            V3 = obj.source.getVertex(i-1,j-1);
            
            [u,v] = obj.getSmoothWeight(V1,V2,V3);
            coordv1 = i*obj.width+j;
            coordv2 = (i-1)*obj.width+j;
            coordv3 = (i-1)*obj.width+j-1;
            
            coordv1 = coordv1 + 1;
            coordv2 = coordv2 + 1;
            coordv3 = coordv3 + 1;
            
            
            %             obj.SmoothCons.add_row(rowCount,obj.x_index(coordv2),(1.0-u)*weight);
            %             obj.SmoothCons.add_row(rowCount,obj.x_index(coordv3),u*weight);
            %             obj.SmoothCons.add_row(rowCount,obj.y_index(coordv3),v*weight);
            %             obj.SmoothCons.add_row(rowCount,obj.y_index(coordv2),(-1.0*v)*weight);
            %             obj.SmoothCons.add_row(rowCount,obj.x_index(coordv1),(-1.0)*weight);
            %             rowCount = rowCount + 1;
            
            
            %V1.x =  V2.x + u * (V3.x - V2.x) + v * (V3.y - V2.y);
            obj.SmoothConstraints(obj.SCc,1) = obj.rowCount; obj.SmoothConstraints(obj.SCc,2) = obj.x_index(coordv2); obj.SmoothConstraints(obj.SCc,3) = (1.0-u)*weight;  obj.SCc = obj.SCc + 1;
            obj.SmoothConstraints(obj.SCc,1) = obj.rowCount; obj.SmoothConstraints(obj.SCc,2) = obj.x_index(coordv3); obj.SmoothConstraints(obj.SCc,3) = u*weight;        obj.SCc = obj.SCc + 1;
            obj.SmoothConstraints(obj.SCc,1) = obj.rowCount; obj.SmoothConstraints(obj.SCc,2) = obj.y_index(coordv3); obj.SmoothConstraints(obj.SCc,3) = v*weight;        obj.SCc = obj.SCc + 1;
            obj.SmoothConstraints(obj.SCc,1) = obj.rowCount; obj.SmoothConstraints(obj.SCc,2) = obj.y_index(coordv2); obj.SmoothConstraints(obj.SCc,3) = (-1.0*v)*weight; obj.SCc = obj.SCc + 1;
            obj.SmoothConstraints(obj.SCc,1) = obj.rowCount; obj.SmoothConstraints(obj.SCc,2) = obj.x_index(coordv1); obj.SmoothConstraints(obj.SCc,3) = (-1.0)*weight;   obj.SCc = obj.SCc + 1;
            obj.rowCount = obj.rowCount+1;
            
            
            
            
            %V1.y = V2.y + u * (V3.y - V2.y) + v * (V2.x - V3.x);
            %             obj.SmoothCons.add_row(rowCount,obj.y_index(coordv2),(1.0-u)*weight);
            %             obj.SmoothCons.add_row(rowCount,obj.y_index(coordv3),u*weight);
            %             obj.SmoothCons.add_row(rowCount,obj.x_index(coordv2),v*weight);
            %             obj.SmoothCons.add_row(rowCount,obj.x_index(coordv3),(-1.0*v)*weight);
            %             obj.SmoothCons.add_row(rowCount,obj.y_index(coordv1),(-1.0)*weight);
            %             rowCount = rowCount + 1;
            
            
            obj.SmoothConstraints(obj.SCc,1) = obj.rowCount; obj.SmoothConstraints(obj.SCc,2) = obj.y_index(coordv2); obj.SmoothConstraints(obj.SCc,3) = (1.0-u)*weight;    obj.SCc = obj.SCc + 1;
            obj.SmoothConstraints(obj.SCc,1) = obj.rowCount; obj.SmoothConstraints(obj.SCc,2) = obj.y_index(coordv3); obj.SmoothConstraints(obj.SCc,3) = u*weight;          obj.SCc = obj.SCc + 1;
            obj.SmoothConstraints(obj.SCc,1) = obj.rowCount; obj.SmoothConstraints(obj.SCc,2) = obj.x_index(coordv2); obj.SmoothConstraints(obj.SCc,3) = v*weight;          obj.SCc = obj.SCc + 1;
            obj.SmoothConstraints(obj.SCc,1) = obj.rowCount; obj.SmoothConstraints(obj.SCc,2) = obj.x_index(coordv3); obj.SmoothConstraints(obj.SCc,3) = (-1.0*v)*weight;   obj.SCc = obj.SCc + 1;
            obj.SmoothConstraints(obj.SCc,1) = obj.rowCount; obj.SmoothConstraints(obj.SCc,2) = obj.y_index(coordv1); obj.SmoothConstraints(obj.SCc,3) = (-1.0)*weight;     obj.SCc = obj.SCc + 1;
            obj.rowCount = obj.rowCount+1;
            
            
        end
        function addCoefficient_3(obj,i,j,weight)
            V1 = obj.source.getVertex(i,j);
            V2 = obj.source.getVertex(i-1,j);
            V3 = obj.source.getVertex(i-1,j+1);
            
            [u,v] = obj.getSmoothWeight(V1,V2,V3);
            coordv1 = i*obj.width+j;
            coordv2 = (i-1)*obj.width+j;
            coordv3 = (i-1)*obj.width+j+1;
            
            coordv1 = coordv1 + 1;
            coordv2 = coordv2 + 1;
            coordv3 = coordv3 + 1;
            
            
            %             obj.SmoothCons.add_row(rowCount,obj.x_index(coordv2),(1.0-u)*weight);
            %             obj.SmoothCons.add_row(rowCount,obj.x_index(coordv3),u*weight);
            %             obj.SmoothCons.add_row(rowCount,obj.y_index(coordv2),v*weight);
            %             obj.SmoothCons.add_row(rowCount,obj.y_index(coordv3),(-1.0*v)*weight);
            %             obj.SmoothCons.add_row(rowCount,obj.x_index(coordv1),(-1.0)*weight);
            %             rowCount = rowCount + 1;
            
            %V1.x=V2.x + u * (V3.x - V2.x) + v * (V2.y - V3.y));
            obj.SmoothConstraints(obj.SCc,1) = obj.rowCount; obj.SmoothConstraints(obj.SCc,2) = obj.x_index(coordv2); obj.SmoothConstraints(obj.SCc,3) = (1.0-u)*weight;  obj.SCc = obj.SCc + 1;
            obj.SmoothConstraints(obj.SCc,1) = obj.rowCount; obj.SmoothConstraints(obj.SCc,2) = obj.x_index(coordv3); obj.SmoothConstraints(obj.SCc,3) = u*weight;        obj.SCc = obj.SCc + 1;
            obj.SmoothConstraints(obj.SCc,1) = obj.rowCount; obj.SmoothConstraints(obj.SCc,2) = obj.y_index(coordv2); obj.SmoothConstraints(obj.SCc,3) = v*weight;        obj.SCc = obj.SCc + 1;
            obj.SmoothConstraints(obj.SCc,1) = obj.rowCount; obj.SmoothConstraints(obj.SCc,2) = obj.y_index(coordv3); obj.SmoothConstraints(obj.SCc,3) = (-1.0*v)*weight; obj.SCc = obj.SCc + 1;
            obj.SmoothConstraints(obj.SCc,1) = obj.rowCount; obj.SmoothConstraints(obj.SCc,2) = obj.x_index(coordv1); obj.SmoothConstraints(obj.SCc,3) = (-1.0)*weight;   obj.SCc = obj.SCc + 1;
            obj.rowCount = obj.rowCount+1;
            
            
            
            
            %V1.y = V2.y + u * (V3.y - V2.y) + v * (V3.x - V2.x));
            %             obj.SmoothCons.add_row(rowCount,obj.y_index(coordv2),(1.0-u)*weight);
            %             obj.SmoothCons.add_row(rowCount,obj.y_index(coordv3),u*weight);
            %             obj.SmoothCons.add_row(rowCount,obj.x_index(coordv3),v*weight);
            %             obj.SmoothCons.add_row(rowCount,obj.x_index(coordv2),(-1.0*v)*weight);
            %             obj.SmoothCons.add_row(rowCount,obj.y_index(coordv1),(-1.0)*weight);
            %             rowCount = rowCount + 1;
            
            obj.SmoothConstraints(obj.SCc,1) = obj.rowCount; obj.SmoothConstraints(obj.SCc,2) = obj.y_index(coordv2); obj.SmoothConstraints(obj.SCc,3) = (1.0-u)*weight;    obj.SCc = obj.SCc + 1;
            obj.SmoothConstraints(obj.SCc,1) = obj.rowCount; obj.SmoothConstraints(obj.SCc,2) = obj.y_index(coordv3); obj.SmoothConstraints(obj.SCc,3) = u*weight;          obj.SCc = obj.SCc + 1;
            obj.SmoothConstraints(obj.SCc,1) = obj.rowCount; obj.SmoothConstraints(obj.SCc,2) = obj.x_index(coordv3); obj.SmoothConstraints(obj.SCc,3) = v*weight;          obj.SCc = obj.SCc + 1;
            obj.SmoothConstraints(obj.SCc,1) = obj.rowCount; obj.SmoothConstraints(obj.SCc,2) = obj.x_index(coordv2); obj.SmoothConstraints(obj.SCc,3) = (-1.0*v)*weight;   obj.SCc = obj.SCc + 1;
            obj.SmoothConstraints(obj.SCc,1) = obj.rowCount; obj.SmoothConstraints(obj.SCc,2) = obj.y_index(coordv1); obj.SmoothConstraints(obj.SCc,3) = (-1.0)*weight;     obj.SCc = obj.SCc + 1;
            obj.rowCount = obj.rowCount+1;
            
            
        end
        function addCoefficient_4(obj,i,j,weight)
            %          V3(i-1,j+1)
            %          |
            %    ______|
            %   V1     V2
            %
            V1 = obj.source.getVertex(i,j);
            V2 = obj.source.getVertex(i,j+1);
            V3 = obj.source.getVertex(i-1,j+1);
            
            [u,v] = obj.getSmoothWeight(V1,V2,V3);
            coordv1 = i*obj.width+j;
            coordv2 = i*obj.width+j+1;
            coordv3 = (i-1)*obj.width+j+1;
            
            coordv1 = coordv1 + 1;
            coordv2 = coordv2 + 1;
            coordv3 = coordv3 + 1;
            
            %V1.x =  V2.x + u * (V3.x - V2.x) + v * (V3.y - V2.y);
            %             obj.SmoothCons.add_row(rowCount,obj.x_index(coordv2),(1.0-u)*weight);
            %             obj.SmoothCons.add_row(rowCount,obj.x_index(coordv3),u*weight);
            %             obj.SmoothCons.add_row(rowCount,obj.y_index(coordv3),v*weight);
            %             obj.SmoothCons.add_row(rowCount,obj.y_index(coordv2),(-1.0*v)*weight);
            %             obj.SmoothCons.add_row(rowCount,obj.x_index(coordv1),(-1.0)*weight);
            %             rowCount = rowCount + 1;
            
            %V1.x =  V2.x + u * (V3.x - V2.x) + v * (V3.y - V2.y);
            obj.SmoothConstraints(obj.SCc,1) = obj.rowCount; obj.SmoothConstraints(obj.SCc,2) = obj.x_index(coordv2); obj.SmoothConstraints(obj.SCc,3) = (1.0-u)*weight;  obj.SCc = obj.SCc + 1;
            obj.SmoothConstraints(obj.SCc,1) = obj.rowCount; obj.SmoothConstraints(obj.SCc,2) = obj.x_index(coordv3); obj.SmoothConstraints(obj.SCc,3) = u*weight;        obj.SCc = obj.SCc + 1;
            obj.SmoothConstraints(obj.SCc,1) = obj.rowCount; obj.SmoothConstraints(obj.SCc,2) = obj.y_index(coordv3); obj.SmoothConstraints(obj.SCc,3) = v*weight;        obj.SCc = obj.SCc + 1;
            obj.SmoothConstraints(obj.SCc,1) = obj.rowCount; obj.SmoothConstraints(obj.SCc,2) = obj.y_index(coordv2); obj.SmoothConstraints(obj.SCc,3) = (-1.0*v)*weight; obj.SCc = obj.SCc + 1;
            obj.SmoothConstraints(obj.SCc,1) = obj.rowCount; obj.SmoothConstraints(obj.SCc,2) = obj.x_index(coordv1); obj.SmoothConstraints(obj.SCc,3) = (-1.0)*weight;   obj.SCc = obj.SCc + 1;
            obj.rowCount = obj.rowCount+1;
            
            
            
            %V1.y = V2.y + u * (V3.y - V2.y) + v * (V2.x - V3.x);
            %             obj.SmoothCons.add_row(rowCount,obj.y_index(coordv2),(1.0-u)*weight);
            %             obj.SmoothCons.add_row(rowCount,obj.y_index(coordv3),u*weight);
            %             obj.SmoothCons.add_row(rowCount,obj.x_index(coordv2),v*weight);
            %             obj.SmoothCons.add_row(rowCount,obj.x_index(coordv3),(-1.0*v)*weight);
            %             obj.SmoothCons.add_row(rowCount,obj.y_index(coordv1),(-1.0)*weight);
            %             rowCount = rowCount + 1;
            
            obj.SmoothConstraints(obj.SCc,1) = obj.rowCount; obj.SmoothConstraints(obj.SCc,2) = obj.y_index(coordv2); obj.SmoothConstraints(obj.SCc,3) = (1.0-u)*weight;    obj.SCc = obj.SCc + 1;
            obj.SmoothConstraints(obj.SCc,1) = obj.rowCount; obj.SmoothConstraints(obj.SCc,2) = obj.y_index(coordv3); obj.SmoothConstraints(obj.SCc,3) = u*weight;          obj.SCc = obj.SCc + 1;
            obj.SmoothConstraints(obj.SCc,1) = obj.rowCount; obj.SmoothConstraints(obj.SCc,2) = obj.x_index(coordv2); obj.SmoothConstraints(obj.SCc,3) = v*weight;          obj.SCc = obj.SCc + 1;
            obj.SmoothConstraints(obj.SCc,1) = obj.rowCount; obj.SmoothConstraints(obj.SCc,2) = obj.x_index(coordv3); obj.SmoothConstraints(obj.SCc,3) = (-1.0*v)*weight;   obj.SCc = obj.SCc + 1;
            obj.SmoothConstraints(obj.SCc,1) = obj.rowCount; obj.SmoothConstraints(obj.SCc,2) = obj.y_index(coordv1); obj.SmoothConstraints(obj.SCc,3) = (-1.0)*weight;     obj.SCc = obj.SCc + 1;
            obj.rowCount = obj.rowCount+1;
            
        end
        function addCoefficient_5(obj,i,j,weight)
            
            %             V1   V2
            %             _____
            %                 |
            %                 |
            %                 V3
            V1 = obj.source.getVertex(i,j);
            V2 = obj.source.getVertex(i,j+1);
            V3 = obj.source.getVertex(i+1,j+1);
            
            [u,v] = obj.getSmoothWeight(V1,V2,V3);
            coordv1 = i*obj.width+j;
            coordv2 = i*obj.width+j+1;
            coordv3 = (i+1)*obj.width+j+1;
            
            coordv1 = coordv1 + 1;
            coordv2 = coordv2 + 1;
            coordv3 = coordv3 + 1;
            
            %V1.x =  V2.x + u * (V3.x - V2.x) + v * (V2.y - V3.y);
            %             obj.SmoothCons.add_row(rowCount,obj.x_index(coordv2),(1.0-u)*weight);
            %             obj.SmoothCons.add_row(rowCount,obj.x_index(coordv3),u*weight);
            %             obj.SmoothCons.add_row(rowCount,obj.y_index(coordv2),v*weight);
            %             obj.SmoothCons.add_row(rowCount,obj.y_index(coordv3),(-1.0*v)*weight);
            %             obj.SmoothCons.add_row(rowCount,obj.x_index(coordv1),(-1.0)*weight);
            %             rowCount = rowCount + 1;
            
            obj.SmoothConstraints(obj.SCc,1) = obj.rowCount; obj.SmoothConstraints(obj.SCc,2) = obj.x_index(coordv2); obj.SmoothConstraints(obj.SCc,3) = (1.0-u)*weight;  obj.SCc = obj.SCc + 1;
            obj.SmoothConstraints(obj.SCc,1) = obj.rowCount; obj.SmoothConstraints(obj.SCc,2) = obj.x_index(coordv3); obj.SmoothConstraints(obj.SCc,3) = u*weight;        obj.SCc = obj.SCc + 1;
            obj.SmoothConstraints(obj.SCc,1) = obj.rowCount; obj.SmoothConstraints(obj.SCc,2) = obj.y_index(coordv2); obj.SmoothConstraints(obj.SCc,3) = v*weight;        obj.SCc = obj.SCc + 1;
            obj.SmoothConstraints(obj.SCc,1) = obj.rowCount; obj.SmoothConstraints(obj.SCc,2) = obj.y_index(coordv3); obj.SmoothConstraints(obj.SCc,3) = (-1.0*v)*weight; obj.SCc = obj.SCc + 1;
            obj.SmoothConstraints(obj.SCc,1) = obj.rowCount; obj.SmoothConstraints(obj.SCc,2) = obj.x_index(coordv1); obj.SmoothConstraints(obj.SCc,3) = (-1.0)*weight;   obj.SCc = obj.SCc + 1;
            obj.rowCount = obj.rowCount+1;
            
            
            
            
            %V1.y = V2.y + u * (V3.y - V2.y) + v * (V3.x - V2.x);
            %             obj.SmoothCons.add_row(rowCount,obj.y_index(coordv2),(1.0-u)*weight);
            %             obj.SmoothCons.add_row(rowCount,obj.y_index(coordv3),u*weight);
            %             obj.SmoothCons.add_row(rowCount,obj.x_index(coordv3),v*weight);
            %             obj.SmoothCons.add_row(rowCount,obj.x_index(coordv2),(-1.0*v)*weight);
            %             obj.SmoothCons.add_row(rowCount,obj.y_index(coordv1),(-1.0)*weight);
            %             rowCount = rowCount + 1;
            
            
            
            obj.SmoothConstraints(obj.SCc,1) = obj.rowCount; obj.SmoothConstraints(obj.SCc,2) = obj.y_index(coordv2); obj.SmoothConstraints(obj.SCc,3) = (1.0-u)*weight;    obj.SCc = obj.SCc + 1;
            obj.SmoothConstraints(obj.SCc,1) = obj.rowCount; obj.SmoothConstraints(obj.SCc,2) = obj.y_index(coordv3); obj.SmoothConstraints(obj.SCc,3) = u*weight;          obj.SCc = obj.SCc + 1;
            obj.SmoothConstraints(obj.SCc,1) = obj.rowCount; obj.SmoothConstraints(obj.SCc,2) = obj.x_index(coordv3); obj.SmoothConstraints(obj.SCc,3) = v*weight;          obj.SCc = obj.SCc + 1;
            obj.SmoothConstraints(obj.SCc,1) = obj.rowCount; obj.SmoothConstraints(obj.SCc,2) = obj.x_index(coordv2); obj.SmoothConstraints(obj.SCc,3) = (-1.0*v)*weight;   obj.SCc = obj.SCc + 1;
            obj.SmoothConstraints(obj.SCc,1) = obj.rowCount; obj.SmoothConstraints(obj.SCc,2) = obj.y_index(coordv1); obj.SmoothConstraints(obj.SCc,3) = (-1.0)*weight;     obj.SCc = obj.SCc + 1;
            obj.rowCount = obj.rowCount+1;
            
        end
        function addCoefficient_6(obj,i,j,weight)
            %             V1(i,j)
            %             |
            %             |____
            %             V2   V3(i+1,j+1)
            %             (i+1,j)
            V1 = obj.source.getVertex(i,j);
            V2 = obj.source.getVertex(i+1,j);
            V3 = obj.source.getVertex(i+1,j+1);
            
            [u,v] = obj.getSmoothWeight(V1,V2,V3);
            coordv1 = i*obj.width+j;
            coordv2 = (i+1)*obj.width+j;
            coordv3 = (i+1)*obj.width+j+1;
            
            coordv1 = coordv1 + 1;
            coordv2 = coordv2 + 1;
            coordv3 = coordv3 + 1;
            
            %V1.x =  V2.x + u * (V3.x - V2.x) + v * (V3.y - V2.y);
            %             obj.SmoothCons.add_row(rowCount,obj.x_index(coordv2),(1.0-u)*weight);
            %             obj.SmoothCons.add_row(rowCount,obj.x_index(coordv3),u*weight);
            %             obj.SmoothCons.add_row(rowCount,obj.y_index(coordv3),v*weight);
            %             obj.SmoothCons.add_row(rowCount,obj.y_index(coordv2),(-1.0*v)*weight);
            %             obj.SmoothCons.add_row(rowCount,obj.x_index(coordv1),(-1.0)*weight);
            %             rowCount = rowCount + 1;
            
            
            obj.SmoothConstraints(obj.SCc,1) = obj.rowCount; obj.SmoothConstraints(obj.SCc,2) = obj.x_index(coordv2); obj.SmoothConstraints(obj.SCc,3) = (1.0-u)*weight;  obj.SCc = obj.SCc + 1;
            obj.SmoothConstraints(obj.SCc,1) = obj.rowCount; obj.SmoothConstraints(obj.SCc,2) = obj.x_index(coordv3); obj.SmoothConstraints(obj.SCc,3) = u*weight;        obj.SCc = obj.SCc + 1;
            obj.SmoothConstraints(obj.SCc,1) = obj.rowCount; obj.SmoothConstraints(obj.SCc,2) = obj.y_index(coordv3); obj.SmoothConstraints(obj.SCc,3) = v*weight;        obj.SCc = obj.SCc + 1;
            obj.SmoothConstraints(obj.SCc,1) = obj.rowCount; obj.SmoothConstraints(obj.SCc,2) = obj.y_index(coordv2); obj.SmoothConstraints(obj.SCc,3) = (-1.0*v)*weight; obj.SCc = obj.SCc + 1;
            obj.SmoothConstraints(obj.SCc,1) = obj.rowCount; obj.SmoothConstraints(obj.SCc,2) = obj.x_index(coordv1); obj.SmoothConstraints(obj.SCc,3) = (-1.0)*weight;   obj.SCc = obj.SCc + 1;
            obj.rowCount = obj.rowCount+1;
            
            
            
            %V1.y = V2.y + u * (V3.y - V2.y) + v * (V2.x - V3.x);
            %             obj.SmoothCons.add_row(rowCount,obj.y_index(coordv2),(1.0-u)*weight);
            %             obj.SmoothCons.add_row(rowCount,obj.y_index(coordv3),u*weight);
            %             obj.SmoothCons.add_row(rowCount,obj.x_index(coordv2),v*weight);
            %             obj.SmoothCons.add_row(rowCount,obj.x_index(coordv3),(-1.0*v)*weight);
            %             obj.SmoothCons.add_row(rowCount,obj.y_index(coordv1),(-1.0)*weight);
            %             rowCount = rowCount + 1;
            
            
            obj.SmoothConstraints(obj.SCc,1) = obj.rowCount; obj.SmoothConstraints(obj.SCc,2) = obj.y_index(coordv2); obj.SmoothConstraints(obj.SCc,3) = (1.0-u)*weight;    obj.SCc = obj.SCc + 1;
            obj.SmoothConstraints(obj.SCc,1) = obj.rowCount; obj.SmoothConstraints(obj.SCc,2) = obj.y_index(coordv3); obj.SmoothConstraints(obj.SCc,3) = u*weight;          obj.SCc = obj.SCc + 1;
            obj.SmoothConstraints(obj.SCc,1) = obj.rowCount; obj.SmoothConstraints(obj.SCc,2) = obj.x_index(coordv2); obj.SmoothConstraints(obj.SCc,3) = v*weight;          obj.SCc = obj.SCc + 1;
            obj.SmoothConstraints(obj.SCc,1) = obj.rowCount; obj.SmoothConstraints(obj.SCc,2) = obj.x_index(coordv3); obj.SmoothConstraints(obj.SCc,3) = (-1.0*v)*weight;   obj.SCc = obj.SCc + 1;
            obj.SmoothConstraints(obj.SCc,1) = obj.rowCount; obj.SmoothConstraints(obj.SCc,2) = obj.y_index(coordv1); obj.SmoothConstraints(obj.SCc,3) = (-1.0)*weight;     obj.SCc = obj.SCc + 1;
            obj.rowCount = obj.rowCount+1;
            
        end
        function addCoefficient_7(obj,i,j,weight)
            %                   V1(i,j)
            %                    |
            %                    |
            %             _______|
            %             V3      V2
            %             (i+1,j-1)
            V1 = obj.source.getVertex(i,j);
            V2 = obj.source.getVertex(i+1,j);
            V3 = obj.source.getVertex(i+1,j-1);
            
            [u,v] = obj.getSmoothWeight(V1,V2,V3);
            coordv1 = i*obj.width+j;
            coordv2 = (i+1)*obj.width+j;
            coordv3 = (i+1)*obj.width+j-1;
            
            coordv1 = coordv1 + 1;
            coordv2 = coordv2 + 1;
            coordv3 = coordv3 + 1;
            
            %V1.x =  V2.x + u * (V3.x - V2.x) + v * (V2.y - V3.y);
            %             obj.SmoothCons.add_row(rowCount,obj.x_index(coordv2),(1.0-u)*weight);
            %             obj.SmoothCons.add_row(rowCount,obj.x_index(coordv3),u*weight);
            %             obj.SmoothCons.add_row(rowCount,obj.y_index(coordv2),v*weight);
            %             obj.SmoothCons.add_row(rowCount,obj.y_index(coordv3),(-1.0*v)*weight);
            %             obj.SmoothCons.add_row(rowCount,obj.x_index(coordv1),(-1.0)*weight);
            %             rowCount = rowCount + 1;
            
            obj.SmoothConstraints(obj.SCc,1) = obj.rowCount; obj.SmoothConstraints(obj.SCc,2) = obj.x_index(coordv2); obj.SmoothConstraints(obj.SCc,3) = (1.0-u)*weight;  obj.SCc = obj.SCc + 1;
            obj.SmoothConstraints(obj.SCc,1) = obj.rowCount; obj.SmoothConstraints(obj.SCc,2) = obj.x_index(coordv3); obj.SmoothConstraints(obj.SCc,3) = u*weight;        obj.SCc = obj.SCc + 1;
            obj.SmoothConstraints(obj.SCc,1) = obj.rowCount; obj.SmoothConstraints(obj.SCc,2) = obj.y_index(coordv2); obj.SmoothConstraints(obj.SCc,3) = v*weight;        obj.SCc = obj.SCc + 1;
            obj.SmoothConstraints(obj.SCc,1) = obj.rowCount; obj.SmoothConstraints(obj.SCc,2) = obj.y_index(coordv3); obj.SmoothConstraints(obj.SCc,3) = (-1.0*v)*weight; obj.SCc = obj.SCc + 1;
            obj.SmoothConstraints(obj.SCc,1) = obj.rowCount; obj.SmoothConstraints(obj.SCc,2) = obj.x_index(coordv1); obj.SmoothConstraints(obj.SCc,3) = (-1.0)*weight;   obj.SCc = obj.SCc + 1;
            obj.rowCount = obj.rowCount+1;
            
            
            
            %V1.y = V2.y + u * (V3.y - V2.y) + v * (V3.x - V2.x);
            %             obj.SmoothCons.add_row(rowCount,obj.y_index(coordv2),(1.0-u)*weight);
            %             obj.SmoothCons.add_row(rowCount,obj.y_index(coordv3),u*weight);
            %             obj.SmoothCons.add_row(rowCount,obj.x_index(coordv3),v*weight);
            %             obj.SmoothCons.add_row(rowCount,obj.x_index(coordv2),(-1.0*v)*weight);
            %             obj.SmoothCons.add_row(rowCount,obj.y_index(coordv1),(-1.0)*weight);
            %             rowCount = rowCount + 1;
            
            
            obj.SmoothConstraints(obj.SCc,1) = obj.rowCount; obj.SmoothConstraints(obj.SCc,2) = obj.y_index(coordv2); obj.SmoothConstraints(obj.SCc,3) = (1.0-u)*weight;    obj.SCc = obj.SCc + 1;
            obj.SmoothConstraints(obj.SCc,1) = obj.rowCount; obj.SmoothConstraints(obj.SCc,2) = obj.y_index(coordv3); obj.SmoothConstraints(obj.SCc,3) = u*weight;          obj.SCc = obj.SCc + 1;
            obj.SmoothConstraints(obj.SCc,1) = obj.rowCount; obj.SmoothConstraints(obj.SCc,2) = obj.x_index(coordv3); obj.SmoothConstraints(obj.SCc,3) = v*weight;          obj.SCc = obj.SCc + 1;
            obj.SmoothConstraints(obj.SCc,1) = obj.rowCount; obj.SmoothConstraints(obj.SCc,2) = obj.x_index(coordv2); obj.SmoothConstraints(obj.SCc,3) = (-1.0*v)*weight;   obj.SCc = obj.SCc + 1;
            obj.SmoothConstraints(obj.SCc,1) = obj.rowCount; obj.SmoothConstraints(obj.SCc,2) = obj.y_index(coordv1); obj.SmoothConstraints(obj.SCc,3) = (-1.0)*weight;     obj.SCc = obj.SCc + 1;
            obj.rowCount = obj.rowCount+1;
            
        end
        function addCoefficient_8(obj,i,j,weight)
            
            %             V2        V1(i,j)
            %             _________
            %             |
            %             |
            %             |
            %             V3(i+1,j-1)
            
            V1 = obj.source.getVertex(i,j);
            V2 = obj.source.getVertex(i,j-1);
            V3 = obj.source.getVertex(i+1,j-1);
            
            [u,v]=obj.getSmoothWeight(V1,V2,V3);
            coordv1 = i*obj.width+j;
            coordv2 = i*obj.width+j-1;
            coordv3 = (i+1)*obj.width+j-1;
            
            coordv1 = coordv1 + 1;
            coordv2 = coordv2 + 1;
            coordv3 = coordv3 + 1;
            
            %V1.x =  V2.x + u * (V3.x - V2.x) + v * (V3.y - V2.y);
            %             obj.SmoothCons.add_row(rowCount,obj.x_index(coordv2),(1.0-u)*weight);
            %             obj.SmoothCons.add_row(rowCount,obj.x_index(coordv3),u*weight);
            %             obj.SmoothCons.add_row(rowCount,obj.y_index(coordv3),v*weight);
            %             obj.SmoothCons.add_row(rowCount,obj.y_index(coordv2),(-1.0*v)*weight);
            %             obj.SmoothCons.add_row(rowCount,obj.x_index(coordv1),(-1.0)*weight);
            %             rowCount = rowCount + 1;
            
            obj.SmoothConstraints(obj.SCc,1) = obj.rowCount; obj.SmoothConstraints(obj.SCc,2) = obj.x_index(coordv2); obj.SmoothConstraints(obj.SCc,3) = (1.0-u)*weight;  obj.SCc = obj.SCc + 1;
            obj.SmoothConstraints(obj.SCc,1) = obj.rowCount; obj.SmoothConstraints(obj.SCc,2) = obj.x_index(coordv3); obj.SmoothConstraints(obj.SCc,3) = u*weight;        obj.SCc = obj.SCc + 1;
            obj.SmoothConstraints(obj.SCc,1) = obj.rowCount; obj.SmoothConstraints(obj.SCc,2) = obj.y_index(coordv3); obj.SmoothConstraints(obj.SCc,3) = v*weight;        obj.SCc = obj.SCc + 1;
            obj.SmoothConstraints(obj.SCc,1) = obj.rowCount; obj.SmoothConstraints(obj.SCc,2) = obj.y_index(coordv2); obj.SmoothConstraints(obj.SCc,3) = (-1.0*v)*weight; obj.SCc = obj.SCc + 1;
            obj.SmoothConstraints(obj.SCc,1) = obj.rowCount; obj.SmoothConstraints(obj.SCc,2) = obj.x_index(coordv1); obj.SmoothConstraints(obj.SCc,3) = (-1.0)*weight;   obj.SCc = obj.SCc + 1;
            obj.rowCount = obj.rowCount+1;
            
            
            %V1.y = V2.y + u * (V3.y - V2.y) + v * (V3.x - V2.x);
            %             obj.SmoothCons.add_row(rowCount,obj.y_index(coordv2),(1.0-u)*weight);
            %             obj.SmoothCons.add_row(rowCount,obj.y_index(coordv3),u*weight);
            %             obj.SmoothCons.add_row(rowCount,obj.x_index(coordv2),v*weight);
            %             obj.SmoothCons.add_row(rowCount,obj.x_index(coordv3),(-1.0*v)*weight);
            %             obj.SmoothCons.add_row(rowCount,obj.y_index(coordv1),(-1.0)*weight);
            %             rowCount = rowCount + 1;
            
            obj.SmoothConstraints(obj.SCc,1) = obj.rowCount; obj.SmoothConstraints(obj.SCc,2) = obj.y_index(coordv2); obj.SmoothConstraints(obj.SCc,3) = (1.0-u)*weight;    obj.SCc = obj.SCc + 1;
            obj.SmoothConstraints(obj.SCc,1) = obj.rowCount; obj.SmoothConstraints(obj.SCc,2) = obj.y_index(coordv3); obj.SmoothConstraints(obj.SCc,3) = u*weight;          obj.SCc = obj.SCc + 1;
            obj.SmoothConstraints(obj.SCc,1) = obj.rowCount; obj.SmoothConstraints(obj.SCc,2) = obj.x_index(coordv2); obj.SmoothConstraints(obj.SCc,3) = v*weight;          obj.SCc = obj.SCc + 1;
            obj.SmoothConstraints(obj.SCc,1) = obj.rowCount; obj.SmoothConstraints(obj.SCc,2) = obj.x_index(coordv3); obj.SmoothConstraints(obj.SCc,3) = (-1.0*v)*weight;   obj.SCc = obj.SCc + 1;
            obj.SmoothConstraints(obj.SCc,1) = obj.rowCount; obj.SmoothConstraints(obj.SCc,2) = obj.y_index(coordv1); obj.SmoothConstraints(obj.SCc,3) = (-1.0)*weight;     obj.SCc = obj.SCc + 1;
            obj.rowCount = obj.rowCount+1;
            
        end       
        function alpha = getAlpha(obj,i,j,weight)            
            count = (obj.histogram(i,j) - obj.threshold)/obj.threshold;
            alpha = weight - count;
            alpha = min(5, alpha);
            alpha = max(0.2, alpha);
            if obj.ADAPTIVE_WEIGHT == 0
                alpha = weight;
            end
        end
    end
    
    
    methods (Access = private) %functions of warping
        function quadWarp(obj,im,q1,q2)
            
            minx = q2.getMinX();
            maxx = q2.getMaxX();
            miny = q2.getMinY();
            maxy = q2.getMaxY();
            
            source = zeros(4,2);
            target = zeros(4,2);
            
            source(1,1) = q2.V00.x;source(1,2) = q2.V00.y;
            source(2,1) = q2.V01.x;source(2,2) = q2.V01.y;
            source(3,1) = q2.V10.x;source(3,2) = q2.V10.y;
            source(4,1) = q2.V11.x;source(4,2) = q2.V11.y;
            
            target(1,1) = q1.V00.x;target(1,2) = q1.V00.y;
            target(2,1) = q1.V01.x;target(2,2) = q1.V01.y;
            target(3,1) = q1.V10.x;target(3,2) = q1.V10.y;
            target(4,1) = q1.V11.x;target(4,2) = q1.V11.y;
            
            H = homography_4pts(source',target');
            H = H./H(3,3);
            
            %qd = Quad(q2.V00,q2.V01,q2.V10,q2.V11);
            obj.warpIm = myWarp(minx,maxx,miny,maxy,im,obj.warpIm,H,obj.gap);
            obj.warpIm = uint8(obj.warpIm);
            
        end
    end
    
    
end

