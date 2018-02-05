classdef VideoStitch2 < handle
    %VIDEOSTITCH Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % inputs
        seqA;
        seqB;
        
        % Dimensions
        nFrames;
        
        videoHeight;
        videoWidth;
        
        meshSize;
        
        quadHeight;
        quadWidth;
        
        % correspondence
        maxppf;
        ppf; % points per frame
        CP; % Control Points
        CPgrid;
        validCP; 
        nCP;
        nCPgrid;
        CPthreshold;
        
        % optimization parameters
        span;
        smoothness;
        cropping;
        stitchness;
        stableW;
        gap;
        
        % Optimization data
        Pa;
        Pb;
        Pa_inv;
        Pb_inv;
        Ca;
        Cb;
        Ca_inv;
        Cb_inv;
        H;
        
        Offset;

        CaCP;
        CbCP;
        
        w_a;
        w_b;
        gamma_a;
        gamma_b;
      
       
    end
    
    methods
        function obj = VideoStitch2(seqA, seqB, pathA, pathB, ControlPoints, ppf, smoothness, cropping, stitchness)
            obj.seqA = seqA;
            obj.seqB = seqB;
            obj.Ca = pathA;
            obj.Cb = pathB;
            obj.Pa = obj.Ca;
            obj.Pb = obj.Cb;
            obj.CP = ControlPoints;
            obj.smoothness = smoothness;
            obj.cropping = cropping;
            obj.stitchness = stitchness;
            obj.stableW = 20;
            obj.span = 30;
            obj.H = eye(3);
            obj.Offset = eye(3);
		
			fileList = dir(seqA);
			fileList = fileList(3:length(fileList));
			obj.nFrames = length(fileList);
			if obj.nFrames < 2
				error('Wrong inputs directory') ;
			end
			fileList = dir(seqB);
			fileList = fileList(3:length(fileList));
			if obj.nFrames ~= length(fileList)
				error('Input length  doesn"t match') ;
			end
			frame = imread([seqB fileList(1).name]);
		
            [obj.videoHeight, obj.videoWidth, ~] = size(frame);
            [~, obj.meshSize, ~, ~, ~] = size(pathA);
            [~, obj.meshSize, ~, ~, ~] = size(pathB);
            [obj.nFrames, obj.maxppf, ~] = size(obj.CP);
            
            % convert CP to CPgrid
            obj.CPgrid = zeros(obj.nFrames, obj.meshSize, obj.meshSize, obj.maxppf / 10, 2);
            obj.nCPgrid = zeros(obj.nFrames, obj.meshSize, obj.meshSize, 2);
                       
            obj.ppf = ppf;
            obj.validCP = zeros(obj.nFrames, obj.maxppf);
            obj.nCP = sum(ppf);
            obj.CPthreshold = zeros(obj.nFrames, 1);
            obj.quadHeight = obj.videoHeight / obj.meshSize;
            obj.quadWidth = obj.videoWidth / obj.meshSize;
            
            for frameIndex = 1:obj.nFrames
                pa = squeeze(obj.CP(frameIndex, 1:obj.ppf(frameIndex), 1:2));
                pb = squeeze(obj.CP(frameIndex, 1:obj.ppf(frameIndex), 3:4));
                
                pacol = floor(pa(:, 1) / obj.quadWidth) + 1;
                parow = floor(pa(:, 2) / obj.quadHeight) + 1;
                pbcol = floor(pb(:, 1) / obj.quadWidth) + 1;
                pbrow = floor(pb(:, 2) / obj.quadHeight) + 1;
                
                for CPindex = 1:obj.ppf(frameIndex)
                    ra = parow(CPindex); ca = pacol(CPindex);
                    rb = pbrow(CPindex); cb = pbcol(CPindex);
                    
                    obj.nCPgrid(frameIndex, ra, ca, 1) = obj.nCPgrid(frameIndex, ra, ca, 1) + 1;
                    obj.CPgrid(frameIndex, ra, ca, obj.nCPgrid(frameIndex, ra, ca, 1), 1) = CPindex;
                    
                    obj.nCPgrid(frameIndex, rb, cb, 2) = obj.nCPgrid(frameIndex, rb, cb, 2) + 1;
                    obj.CPgrid(frameIndex, rb, cb, obj.nCPgrid(frameIndex, rb, cb, 2), 2) = CPindex;
                end            
            end            
            %compute C_inv
            obj.Ca_inv = zeros(size(obj.Ca));
            obj.Cb_inv = zeros(size(obj.Cb));
            
            
            obj.CaCP = zeros(obj.nFrames, obj.maxppf, 3);
            obj.CbCP = obj.CaCP;
            
            obj.w_a = zeros(obj.nFrames, obj.nFrames, obj.meshSize, obj.meshSize);
            obj.w_b = zeros(obj.nFrames, obj.nFrames, obj.meshSize, obj.meshSize);
            
            obj.gamma_a = zeros(obj.nFrames, obj.meshSize, obj.meshSize);
            obj.gamma_b = zeros(obj.nFrames, obj.meshSize, obj.meshSize);
        end
        
        function [Pa_inv, Pb_inv] = inverse(obj, Pa, Pb)
            Pa_inv = zeros(size(Pa));
            Pb_inv = zeros(size(Pb));
            for frameIndex = 1:obj.nFrames
                for i = 1:obj.meshSize
                    for j = 1:obj.meshSize
                        Pa_inv(frameIndex, i, j, :, :) = squeeze(Pa(frameIndex, i, j, :, :))^(-1);                         
                        Pb_inv(frameIndex, i, j, :, :) = squeeze(Pb(frameIndex, i, j, :, :))^(-1);   
                    end
                end
            end
        end
        
        function [PaCorner, PbCorner] = PxCorner(obj, Pa_inv, Pb_inv)
            PaCorner = ones(obj.nFrames, obj.meshSize, obj.meshSize, 4, 3);
            PbCorner = ones(obj.nFrames, obj.meshSize, obj.meshSize, 4, 3);
            for frameIndex = 1:obj.nFrames
                for i = 1:obj.meshSize
                    for j = 1:obj.meshSize
                        [PaCorner(frameIndex, i, j, 1, 1), PaCorner(frameIndex, i, j, 1, 2)]...
                            = obj.transform([(j - 1) * obj.quadWidth + 1 (i - 1) * obj.quadHeight + 1], squeeze(Pa_inv(frameIndex, i, j, :, :))); 
                        [PaCorner(frameIndex, i, j, 2, 1), PaCorner(frameIndex, i, j, 2, 2)]...
                            = obj.transform([(j) * obj.quadWidth (i - 1) * obj.quadHeight + 1], squeeze(Pa_inv(frameIndex, i, j, :, :)));
                        [PaCorner(frameIndex, i, j, 3, 1), PaCorner(frameIndex, i, j, 3, 2)]...
                            = obj.transform([(j - 1) * obj.quadWidth + 1 (i) * obj.quadHeight], squeeze(Pa_inv(frameIndex, i, j, :, :)));
                        [PaCorner(frameIndex, i, j, 4, 1), PaCorner(frameIndex, i, j, 4, 2)]...
                            = obj.transform([(j) * obj.quadWidth (i) * obj.quadHeight], squeeze(Pa_inv(frameIndex, i, j, :, :)));
                        [PbCorner(frameIndex, i, j, 1, 1), PbCorner(frameIndex, i, j, 1, 2)]...
                            = obj.transform([(j - 1) * obj.quadWidth + 1 (i - 1) * obj.quadHeight + 1], squeeze(Pb_inv(frameIndex, i, j, :, :))); 
                        [PbCorner(frameIndex, i, j, 2, 1), PbCorner(frameIndex, i, j, 2, 2)]...
                            = obj.transform([(j) * obj.quadWidth (i - 1) * obj.quadHeight + 1], squeeze(Pb_inv(frameIndex, i, j, :, :)));
                        [PbCorner(frameIndex, i, j, 3, 1), PbCorner(frameIndex, i, j, 3, 2)]...
                            = obj.transform([(j - 1) * obj.quadWidth + 1 (i) * obj.quadHeight], squeeze(Pb_inv(frameIndex, i, j, :, :)));
                        [PbCorner(frameIndex, i, j, 4, 1), PbCorner(frameIndex, i, j, 4, 2)]...
                            = obj.transform([(j) * obj.quadWidth (i) * obj.quadHeight], squeeze(Pb_inv(frameIndex, i, j, :, :)));
                    end
                end
            end
        end
        
        function init(obj)
            obj.Pa = obj.Ca;
            obj.Pb = obj.Cb;
            for frameIndex = 1:obj.nFrames
                for k = 1:obj.ppf(frameIndex)
                    obj.validCP(frameIndex, k) = 1;
                end
            end
            obj.calcOmega();
            obj.H = eye(3);
            obj.nCP = sum(obj.ppf);
            fprintf('Number of Valid Control Points :%5d\n', obj.nCP);
        end
        
        function updateOffset(obj)
            B = zeros(obj.nFrames, 3, 3);
            ms = obj.meshSize;
            for frameIndex = 1:obj.nFrames
                for row = 1:ms
                    for col = 1:ms 
                        B1 = squeeze(obj.Pa(frameIndex, row, col, :, :)) / squeeze(obj.Ca(frameIndex, row, col, :, :));
                        B(frameIndex, :, :) = squeeze(B(frameIndex, :, :)) + B1 ./ B1(3,3); 
                    end
                end
            end
            nf = obj.nFrames;
            offset = squeeze(sum(B(:, :, :)) / nf / (ms * ms));
%             offset = eye(3);
            t = floor(obj.nFrames * 0.5);
            [ax00, ay00] = obj.transform([1, 1], offset \ squeeze(obj.Pa(t, obj.meshSize, 1, :, :)) / squeeze(obj.Ca(t, obj.meshSize, 1, :, :)));
            [ax01, ay01] = obj.transform([1, obj.videoHeight], offset \ squeeze(obj.Pa(t, 1, 1, :, :)) / squeeze(obj.Ca(t, 1, 1, :, :)));
            [ax10, ay10] = obj.transform([obj.videoWidth, 1], offset \ squeeze(obj.Pa(t, obj.meshSize, obj.meshSize, :, :)) / squeeze(obj.Ca(t, obj.meshSize, obj.meshSize, :, :)));
            [ax11, ay11] = obj.transform([obj.videoWidth, obj.videoHeight], offset \ squeeze(obj.Pa(t, obj.meshSize, 1, :, :)) / squeeze(obj.Ca(t, obj.meshSize, 1, :, :)));
            
            [bx00, by00] = obj.transform([1, 1], offset \ obj.H * squeeze(obj.Pb(t, obj.meshSize, 1, :, :)) / squeeze(obj.Cb(t, obj.meshSize, 1, :, :)));
            [bx01, by01] = obj.transform([1, obj.videoHeight], offset \ obj.H * squeeze(obj.Pb(t, 1, 1, :, :)) / squeeze(obj.Cb(t, 1, 1, :, :)));
            [bx10, by10] = obj.transform([obj.videoWidth, 1], offset \ obj.H * squeeze(obj.Pb(t, obj.meshSize, obj.meshSize, :, :)) / squeeze(obj.Cb(t, obj.meshSize, obj.meshSize, :, :)));
            [bx11, by11] = obj.transform([obj.videoWidth, obj.videoHeight], offset \ obj.H * squeeze(obj.Pb(t, obj.meshSize, 1, :, :)) / squeeze(obj.Cb(t, obj.meshSize, 1, :, :)));
            
            a = [ax00 ay00; ax01 ay01; ax10 ay10; ax11 ay11];
            b = [bx00 by00; bx01 by01; bx10 by10; bx11 by11];
            c = (a + b) * 0.5;
            
            ha = homography_4pts(c', a');
            hb = homography_4pts(b', c');
            obj.Offset = ha / offset;
        end
        
        function value = getCoherenceTerm(obj, ab, P, row, col, frameIndex)
            value = 0;
            if ab == 'a'
                C = obj.Ca;
            else
                C = obj.Cb;                
            end
            if obj.meshSize == 1
                return 
            end
            
            % i=1, j=1
            if row == 1 && col == 1
                value = value + 2 * obj.stableW * squeeze(P(frameIndex, row+1, col, :, :)) / squeeze(C(frameIndex, row+1, col, :, :)) * squeeze(C(frameIndex, row, col, :, :));
                value = value + 2 * obj.stableW * squeeze(P(frameIndex, row, col+1, :, :)) / squeeze(C(frameIndex, row, col+1, :, :)) * squeeze(C(frameIndex, row, col, :, :));
                value = value + 2 * obj.stableW * squeeze(P(frameIndex, row+1, col+1, :, :)) / squeeze(C(frameIndex, row+1, col+1, :, :)) * squeeze(C(frameIndex, row, col, :, :));
            end
            % row=1, col=2:15
            if row == 1 && col < obj.meshSize && col > 1
                value = value + 2 * obj.stableW * squeeze(P(frameIndex, row+1, col, :, :)) / squeeze(C(frameIndex, row+1, col, :, :)) * squeeze(C(frameIndex, row, col, :, :));
                value = value + 2 * obj.stableW * squeeze(P(frameIndex, row, col+1, :, :)) / squeeze(C(frameIndex, row, col+1, :, :)) * squeeze(C(frameIndex, row, col, :, :));
                value = value + 2 * obj.stableW * squeeze(P(frameIndex, row+1, col+1, :, :)) / squeeze(C(frameIndex, row+1, col+1, :, :)) * squeeze(C(frameIndex, row, col, :, :));
                value = value + 2 * obj.stableW * squeeze(P(frameIndex, row+1, col-1, :, :)) / squeeze(C(frameIndex, row+1, col-1, :, :)) * squeeze(C(frameIndex, row, col, :, :));
                value = value + 2 * obj.stableW * squeeze(P(frameIndex, row, col-1, :, :)) / squeeze(C(frameIndex, row, col-1, :, :)) * squeeze(C(frameIndex, row, col, :, :));
            end
            % row=1, col=16
            if row == 1 && col == obj.meshSize
                value = value + 2 * obj.stableW * squeeze(P(frameIndex, row+1, col, :, :)) / squeeze(C(frameIndex, row+1, col, :, :)) * squeeze(C(frameIndex, row, col, :, :));
                value = value + 2 * obj.stableW * squeeze(P(frameIndex, row, col-1, :, :)) / squeeze(C(frameIndex, row, col-1, :, :)) * squeeze(C(frameIndex, row, col, :, :));
                value = value + 2 * obj.stableW * squeeze(P(frameIndex, row+1, col-1, :, :)) / squeeze(C(frameIndex, row+1, col-1, :, :)) * squeeze(C(frameIndex, row, col, :, :));
            end
            % row=2:15, col=1
            if col == 1 && row < obj.meshSize && row > 1
                value = value + 2 * obj.stableW * squeeze(P(frameIndex, row+1, col, :, :)) / squeeze(C(frameIndex, row+1, col, :, :)) * squeeze(C(frameIndex, row, col, :, :));
                value = value + 2 * obj.stableW * squeeze(P(frameIndex, row-1, col, :, :)) / squeeze(C(frameIndex, row-1, col, :, :)) * squeeze(C(frameIndex, row, col, :, :));
                value = value + 2 * obj.stableW * squeeze(P(frameIndex, row+1, col+1, :, :)) / squeeze(C(frameIndex, row+1, col+1, :, :)) * squeeze(C(frameIndex, row, col, :, :));
                value = value + 2 * obj.stableW * squeeze(P(frameIndex, row-1, col+1, :, :)) / squeeze(C(frameIndex, row-1, col+1, :, :)) * squeeze(C(frameIndex, row, col, :, :));
                value = value + 2 * obj.stableW * squeeze(P(frameIndex, row, col+1, :, :)) / squeeze(C(frameIndex, row, col+1, :, :)) * squeeze(C(frameIndex, row, col, :, :));
            end
            % row=2:15, col=2:15
            if row > 1 && row < obj.meshSize && col > 1 && col < obj.meshSize
                value = value + 2 * obj.stableW * squeeze(P(frameIndex, row+1, col, :, :)) / squeeze(C(frameIndex, row+1, col, :, :)) * squeeze(C(frameIndex, row, col, :, :));
                value = value + 2 * obj.stableW * squeeze(P(frameIndex, row-1, col, :, :)) / squeeze(C(frameIndex, row-1, col, :, :)) * squeeze(C(frameIndex, row, col, :, :));
                value = value + 2 * obj.stableW * squeeze(P(frameIndex, row+1, col+1, :, :)) / squeeze(C(frameIndex, row+1, col+1, :, :)) * squeeze(C(frameIndex, row, col, :, :));
                value = value + 2 * obj.stableW * squeeze(P(frameIndex, row-1, col+1, :, :)) / squeeze(C(frameIndex, row-1, col+1, :, :)) * squeeze(C(frameIndex, row, col, :, :));
                value = value + 2 * obj.stableW * squeeze(P(frameIndex, row-1, col-1, :, :)) / squeeze(C(frameIndex, row-1, col-1, :, :)) * squeeze(C(frameIndex, row, col, :, :));
                value = value + 2 * obj.stableW * squeeze(P(frameIndex, row+1, col-1, :, :)) / squeeze(C(frameIndex, row+1, col-1, :, :)) * squeeze(C(frameIndex, row, col, :, :));
                value = value + 2 * obj.stableW * squeeze(P(frameIndex, row, col+1, :, :)) / squeeze(C(frameIndex, row, col+1, :, :)) * squeeze(C(frameIndex, row, col, :, :));
                value = value + 2 * obj.stableW * squeeze(P(frameIndex, row, col-1, :, :)) / squeeze(C(frameIndex, row, col-1, :, :)) * squeeze(C(frameIndex, row, col, :, :)); 
            end
            % row=2:15, col=16
            if col == obj.meshSize && row < obj.meshSize && row > 1
                value = value + 2 * obj.stableW * squeeze(P(frameIndex, row+1, col, :, :)) / squeeze(C(frameIndex, row+1, col, :, :)) * squeeze(C(frameIndex, row, col, :, :));
                value = value + 2 * obj.stableW * squeeze(P(frameIndex, row, col-1, :, :)) / squeeze(C(frameIndex, row, col-1, :, :)) * squeeze(C(frameIndex, row, col, :, :));
                value = value + 2 * obj.stableW * squeeze(P(frameIndex, row-1, col, :, :)) / squeeze(C(frameIndex, row-1, col, :, :)) * squeeze(C(frameIndex, row, col, :, :));
                value = value + 2 * obj.stableW * squeeze(P(frameIndex, row+1, col-1, :, :)) / squeeze(C(frameIndex, row+1, col-1, :, :)) * squeeze(C(frameIndex, row, col, :, :));
                value = value + 2 * obj.stableW * squeeze(P(frameIndex, row-1, col-1, :, :)) / squeeze(C(frameIndex, row-1, col-1, :, :)) * squeeze(C(frameIndex, row, col, :, :));
            end
            % row=16, col=1
            if row == obj.meshSize && col == 1
                value = value + 2 * obj.stableW * squeeze(P(frameIndex, row-1, col, :, :)) / squeeze(C(frameIndex, row-1, col, :, :)) * squeeze(C(frameIndex, row, col, :, :));
                value = value + 2 * obj.stableW * squeeze(P(frameIndex, row, col+1, :, :)) / squeeze(C(frameIndex, row, col+1, :, :)) * squeeze(C(frameIndex, row, col, :, :));
                value = value + 2 * obj.stableW * squeeze(P(frameIndex, row-1, col+1, :, :)) / squeeze(C(frameIndex, row-1, col+1, :, :)) * squeeze(C(frameIndex, row, col, :, :));
            end
            % row=16, col=2:15
            if row == obj.meshSize && col < obj.meshSize && col > 1
                value = value + 2 * obj.stableW * squeeze(P(frameIndex, row-1, col, :, :)) / squeeze(C(frameIndex, row-1, col, :, :)) * squeeze(C(frameIndex, row, col, :, :));
                value = value + 2 * obj.stableW * squeeze(P(frameIndex, row, col+1, :, :)) / squeeze(C(frameIndex, row, col+1, :, :)) * squeeze(C(frameIndex, row, col, :, :));
                value = value + 2 * obj.stableW * squeeze(P(frameIndex, row-1, col+1, :, :)) / squeeze(C(frameIndex, row-1, col+1, :, :)) * squeeze(C(frameIndex, row, col, :, :));
                value = value + 2 * obj.stableW * squeeze(P(frameIndex, row-1, col-1, :, :)) / squeeze(C(frameIndex, row-1, col-1, :, :)) * squeeze(C(frameIndex, row, col, :, :));
                value = value + 2 * obj.stableW * squeeze(P(frameIndex, row, col-1, :, :)) / squeeze(C(frameIndex, row, col-1, :, :)) * squeeze(C(frameIndex, row, col, :, :));
            end
            % row=16, col=16
            if row == obj.meshSize && col == obj.meshSize
                value = value + 2 * obj.stableW * squeeze(P(frameIndex, row-1, col, :, :)) / squeeze(C(frameIndex, row-1, col, :, :)) * squeeze(C(frameIndex, row, col, :, :));
                value = value + 2 * obj.stableW * squeeze(P(frameIndex, row, col-1, :, :)) / squeeze(C(frameIndex, row, col-1, :, :)) * squeeze(C(frameIndex, row, col, :, :));
                value = value + 2 * obj.stableW * squeeze(P(frameIndex, row-1, col-1, :, :)) / squeeze(C(frameIndex, row-1, col-1, :, :)) * squeeze(C(frameIndex, row, col, :, :));
            end
            
        end
        
        function optPath(obj, maxIte, secondPhase)
            % zero round
            firstround = maxIte - secondPhase;
            % parameters
            inMaxIte = 10; % number of innner-iteration for the first outer-iteration
            asap_1 = 0.6;
            asap_2 = 0.3;
            boost_non_commmon = 5;
            boost_commmon = 0.1;
            
            for ite = 1:maxIte
                fprintf('\nRound#%2d ', ite);                
                if ite > 1
                    inMaxIte = 3;
                    if ite > firstround
                        asaplambda = asap_2;
                        inMaxIte = 1;
                        obj.stableW = 30; % 20-40 is OK, larger for more regidity.
                    else
                        asaplambda = asap_1;
                    end
                    tic;
                    parfor frameIndex = 1:obj.nFrames
                        CPcount = 0;
                        PAa = zeros(obj.ppf(frameIndex), 2);
                        PBa = zeros(obj.ppf(frameIndex), 2);
                        PAb = zeros(obj.ppf(frameIndex), 2);
                        PBb = zeros(obj.ppf(frameIndex), 2);
                        Ba = zeros(obj.meshSize, obj.meshSize, 3, 3);
                        Bb = Ba;
                        
                        for row = 1:obj.meshSize
                            for col = 1:obj.meshSize
                                Ba(row, col, :, :) = squeeze(obj.Pa(frameIndex, row, col, :, :)) / squeeze(obj.Ca(frameIndex, row, col, :, :));
                                Bb(row, col, :, :) = squeeze(obj.Pb(frameIndex, row, col, :, :)) / squeeze(obj.Cb(frameIndex, row, col, :, :));
                            end
                        end
                        
                        for k = 1:obj.maxppf
                            if obj.validCP(frameIndex, k) == 0
                                break;
                            end
                            if obj.validCP(frameIndex, k) == 1
                                pa = obj.CP(frameIndex, k, 1:2);
                                cola = floor((pa(1) - 0.001) / obj.quadWidth) + 1;
                                rowa = floor((pa(2) - 0.001) / obj.quadHeight) + 1;
                                pb = obj.CP(frameIndex, k, 3:4);
                                colb = floor((pb(1) - 0.001) / obj.quadWidth) + 1;
                                rowb = floor((pb(2) - 0.001) / obj.quadHeight) + 1;
                                CPcount = CPcount + 1;
                                PAa(CPcount, :) = pa;
                                [PBa(CPcount, 1), PBa(CPcount, 2)] = obj.transform(pb, obj.H * squeeze(Bb(rowb, colb, :, :)));    
                                PBb(CPcount, :) = pb;
                                [PAb(CPcount, 1), PAb(CPcount, 2)] = obj.transform(pa, obj.H \ squeeze(Ba(rowa, cola, :, :)));    
                            end
                        end                            
                        PaStitch(frameIndex, :, :, :, :) = NewWarping(PAa, PBa, obj.videoHeight, obj.videoWidth, obj.quadHeight, obj.quadWidth, asaplambda);
                        PbStitch(frameIndex, :, :, :, :) = NewWarping(PBb, PAb, obj.videoHeight, obj.videoWidth, obj.quadHeight, obj.quadWidth, asaplambda);
                    end
                    toc; 
                end
                obj.calcOmega();
                for inIte = 1:inMaxIte
                    fprintf('.');
                    oPa = obj.Pa;
                    oPb = obj.Pb;                        
                    for frameIndex = 1:obj.nFrames
                        stitching = obj.stitchness;
                        % fix the head and tail
                        crop_fix = 0;
                        if ite == 1
                            if frameIndex <= obj.span
                                crop_fix = 1 * (obj.span - frameIndex + 1);
                            end
                            if frameIndex > obj.nFrames - obj.span
                                crop_fix = 1 * ( - obj.nFrames + frameIndex + obj.span);
                            end
                        end
                        
                        for row = 1:obj.meshSize
                            for col = 1:obj.meshSize
                                head = max(frameIndex - obj.span, 1);
                                tail = min(frameIndex + obj.span, obj.nFrames);
                                nn = tail - head + 1;
                                % Pat - Cat
                                value = squeeze(obj.Pa(frameIndex, row, col, :, :)) * (crop_fix + 1);                            
                                % Pat - Par
                                n9 = reshape(oPa(head:tail, row, col, :, :), [nn, 9]);
                                weight = obj.w_a(frameIndex, head:tail, row, col);
                                value = value + 2 * obj.smoothness * reshape(weight * n9, [3, 3]);
                                % Pat - Pat
                                value = value + obj.getCoherenceTerm('a', oPa, row, col, frameIndex);
                                % Pat - Pbt'
                                
                                obj.Pa(frameIndex, row, col, :, :) = value / (crop_fix + obj.gamma_a(frameIndex, row, col));
                                
                                if ite == 1
                                    obj.Pa(frameIndex, row, col, :, :) = value / (crop_fix + obj.gamma_a(frameIndex, row, col));
                                end
                                if ite <= firstround && ite > 1
                                    value = value + squeeze(PaStitch(frameIndex, row, col, :, :)) * squeeze(obj.Ca(frameIndex, row, col, :, :)) * stitching;
                                    obj.Pa(frameIndex, row, col, :, :) = value / (crop_fix + obj.gamma_a(frameIndex, row, col) + stitching);
                                end
                                if ite > firstround && obj.nCPgrid(frameIndex, row, col, 1) > 10
                                    value = value + squeeze(PaStitch(frameIndex, row, col, :, :)) * squeeze(obj.Ca(frameIndex, row, col, :, :)) * stitching * boost_non_commmon;
                                    obj.Pa(frameIndex, row, col, :, :) = value / (crop_fix + obj.gamma_a(frameIndex, row, col) + stitching * boost_non_commmon);
                                end
                                
                                if ite > firstround && obj.nCPgrid(frameIndex, row, col, 1) <= 10
                                    value = value + squeeze(PaStitch(frameIndex, row, col, :, :)) * squeeze(obj.Ca(frameIndex, row, col, :, :)) * stitching * boost_commmon;
                                    obj.Pa(frameIndex, row, col, :, :) = value / (crop_fix + obj.gamma_a(frameIndex, row, col) + stitching * boost_commmon);
                                end

                                % Pb
                                value = squeeze(obj.Pb(frameIndex, row, col, :, :)) * (crop_fix + 1);
                                n9 = reshape(oPb(head:tail, row, col, :, :), [nn, 9]);
                                weight = obj.w_b(frameIndex, head:tail, row, col);
                                value = value + 2 * obj.smoothness * reshape(weight * n9, [3, 3]);

                                value = value + obj.getCoherenceTerm('b', oPb, row, col, frameIndex);
                                if ite == 1
                                    obj.Pb(frameIndex, row, col, :, :) = value / (crop_fix + obj.gamma_b(frameIndex, row, col));
                                end                                
                                if ite <= firstround && ite > 1 
                                    value = value + squeeze(PbStitch(frameIndex, row, col, :, :)) * squeeze(obj.Cb(frameIndex, row, col, :, :)) * stitching;
                                    obj.Pb(frameIndex, row, col, :, :) = value / (crop_fix + obj.gamma_b(frameIndex, row, col) + stitching);
                                end                                
                                if ite > firstround && obj.nCPgrid(frameIndex, row, col, 2) > 10
                                    value = value + squeeze(PbStitch(frameIndex, row, col, :, :)) * squeeze(obj.Cb(frameIndex, row, col, :, :)) * stitching * boost_non_commmon;                                    
                                    obj.Pb(frameIndex, row, col, :, :) = value / (crop_fix + obj.gamma_b(frameIndex, row, col) + stitching * boost_non_commmon);
                                end
                                if ite > firstround && obj.nCPgrid(frameIndex, row, col, 2) <= 10
                                    value = value + squeeze(PbStitch(frameIndex, row, col, :, :)) * squeeze(obj.Cb(frameIndex, row, col, :, :)) * stitching * boost_commmon;                                    
                                    obj.Pb(frameIndex, row, col, :, :) = value / (crop_fix + obj.gamma_b(frameIndex, row, col) + stitching * boost_commmon);
                                end
                                
                            end
                        end
                    end
                end                
                obj.updateH_simple();
                obj.updateCP();
            end
        end
        
        function updateH_simple(obj)
            % compute H from Pa Ca Pb Cb and CP 
            PA = zeros(2, obj.ppf(obj.span+1) + obj.ppf(obj.span+2));
            PB = zeros(2, obj.ppf(obj.span+1) + obj.ppf(obj.span+2));
            CPcount = 0;
            for frameIndex = floor(obj.nFrames * 0.5) : floor(obj.nFrames * 0.5) + 1
                for k = 1:obj.maxppf
                    if obj.validCP(frameIndex, k) == 0
                        break                        
                    end
                    if obj.validCP(frameIndex, k) == 1
                        %xa, ya, xb, yb;
                        pb = obj.CP(frameIndex, k, 3:4);
                        colb = floor((pb(1) - 0.001) / obj.quadWidth) + 1;
                        rowb = floor((pb(2) - 0.001) / obj.quadHeight) + 1;
                        Bb = squeeze(obj.Pb(frameIndex, rowb, colb, :, :)) / squeeze(obj.Cb(frameIndex, rowb, colb, :, :));
                        pa = obj.CP(frameIndex, k, 1:2);
                        cola = floor((pa(1) - 0.001) / obj.quadWidth) + 1;
                        rowa = floor((pa(2) - 0.001) / obj.quadHeight) + 1; 
                        Ba = squeeze(obj.Pa(frameIndex, rowa, cola, :, :)) / squeeze(obj.Ca(frameIndex, rowa, cola, :, :));
                        CPcount = CPcount + 1;
                        [PA(1, CPcount), PA(2, CPcount)] = obj.transform(pa, Ba);
                        [PB(1, CPcount), PB(2, CPcount)] = obj.transform(pb, Bb);
                    end                        
                end
            end
            obj.H = ransacfithomography(PA, PB, 0.001);
            obj.H = obj.H / obj.H(3, 3);
            obj.H = obj.H ^ (-1);
        end
        
        function updateCP(obj)
            % set validCp and ppf according to H
            obj.ppf = zeros(obj.nFrames, 1);
            for frameIndex = 1:obj.nFrames
                d = zeros(obj.maxppf, 1);
                for k = 1:obj.maxppf
                    if obj.validCP(frameIndex, k) == 0
                        break
                    end
                    pb = obj.CP(frameIndex, k, 3:4);
                    colb = floor((pb(1) - 0.001) / obj.quadWidth) + 1;
                    rowb = floor((pb(2) - 0.001) / obj.quadHeight) + 1;
                    pa = obj.CP(frameIndex, k, 1:2);
                    cola = floor((pa(1) - 0.001) / obj.quadWidth) + 1;
                    rowa = floor((pa(2) - 0.001) / obj.quadHeight) + 1;
                    
                    if cola < 1 || cola > obj.meshSize
                        cola = 1;
                    end
                    if colb < 1 || colb > obj.meshSize
                        cola = 1;
                    end
                    if rowa < 1 || rowa > obj.meshSize
                        cola = 1;
                    end
                    if rowb < 1 || rowb > obj.meshSize
                        cola = 1;
                    end
                    
                    [a2, b2] = obj.transform(pb, obj.H * squeeze(obj.Pb(frameIndex, rowb, colb, :, :)) / squeeze(obj.Cb(frameIndex, rowb, colb, :, :)));
                    [a1, b1] = obj.transform(pa, squeeze(obj.Pa(frameIndex, rowa, cola, :, :)) / squeeze(obj.Ca(frameIndex, rowa, cola, :, :)));
                    d(k) = sqrt((a1 - a2)^2 + (b1 - b2)^2);                            
                end
                dsort = sort(d(d ~= 0));
                threshold = mean(dsort(floor(length(dsort) * 2 / 5):floor(length(dsort) * 16/ 20))) * 2 + 10;

                obj.validCP(frameIndex, d <= threshold) = 1;
                obj.validCP(frameIndex, d == 0) = 0;
                obj.validCP(frameIndex, d > threshold) = 0;
                obj.ppf(frameIndex) = sum(obj.validCP(frameIndex, :));
                obj.validCP(frameIndex, d > threshold) = -1;
                if obj.ppf(frameIndex) < 20
                    obj.validCP(frameIndex, obj.validCP(frameIndex, :) ~= 0) = -1; 
                    obj.ppf(frameIndex) = 0;
                end
                obj.CPthreshold(frameIndex) = threshold;
            end
            obj.nCP = sum(obj.ppf);
            fprintf('\nNumber of Valid Control Points :%5d\n', obj.nCP);
            fprintf('\nAverage Alignment Error :%.2f\n', mean(obj.CPthreshold - 10) / 2);
            obj.nCPgrid = zeros(obj.nFrames, obj.meshSize, obj.meshSize, 2);
            for frameIndex = 1:obj.nFrames
                for k = 1:obj.maxppf
                    if obj.validCP(frameIndex, k) ==1
                        pb = obj.CP(frameIndex, k, 3:4);
                        colb = floor((pb(1) - 0.001) / obj.quadWidth) + 1;
                        rowb = floor((pb(2) - 0.001) / obj.quadHeight) + 1;
                        pa = obj.CP(frameIndex, k, 1:2);
                        cola = floor((pa(1) - 0.001) / obj.quadWidth) + 1;
                        rowa = floor((pa(2) - 0.001) / obj.quadHeight) + 1;
                        obj.nCPgrid(frameIndex, rowa, cola, 1) = obj.nCPgrid(frameIndex, rowa, cola, 1) + 1;
                        obj.nCPgrid(frameIndex, rowb, colb, 2) = obj.nCPgrid(frameIndex, rowb, colb, 2) + 1;
                    end                    
                end
            end 
        end
        
        
                
        function render(obj, outPath, gap)
            if ~exist(outPath, 'dir')
                mkdir(outPath);
            end
            obj.gap = gap;
            obj.updateOffset();
            parfor frameIndex = 1 : obj.nFrames %parfor
                disp(['rendering: # ' int2str(frameIndex)]);
                fileListA = dir(obj.seqA);
                fileListA = fileListA(3:length(fileListA));
                fileListB = dir(obj.seqB);
                fileListB = fileListB(3:length(fileListB));
                fileNameA = fileListA(frameIndex).name;
                fileNameB = fileListB(frameIndex).name;
                IA = imread([obj.seqA fileNameA]);
                IB = imread([obj.seqB fileNameB]);              
% ----------------DEBUG_MODE------------------
%                 warpA = zeros(obj.videoHeight + 2 * obj.gap, obj.videoWidth + 2 * obj.gap, 3);
%                 warpB = zeros(obj.videoHeight + 2 * obj.gap, obj.videoWidth + 2 * obj.gap, 3);
%                 for i = 1:obj.meshSize
%                     for j = 1:obj.meshSize
%                         warpA = obj.warp1('a', frameIndex, i, j, IA, warpA); 
%                         warpB = obj.warp1('b', frameIndex, i, j, IB, warpB); 
%                     end
%                 end
                
                warpA = obj.render1('a', IA, frameIndex, obj.Pa, obj.Ca, 0);
                warpB = obj.render1('b', IB, frameIndex, obj.Pb, obj.Cb, 0);
                
                
                imwrite(uint8(warpA), [outPath '/A' int2str(frameIndex) '.png']);
                imwrite(uint8(warpB), [outPath '/B' int2str(frameIndex) '.png']);

                warpA = uint8(warpA); warpB = uint8(warpB);
                warp = obj.imageBlending(warpA, warpB);
                
                imwrite(warp, [outPath '/' int2str(frameIndex) '.png']);
            end
        end
        
        function renderGrid(obj, outPath, gap)
            if ~exist(outPath, 'dir')
                mkdir(outPath);
            end
            obj.gap = gap;
            obj.updateOffset();
            parfor frameIndex = 1:obj.nFrames %parfor
                disp(['rendering: # ' int2str(frameIndex)]);
                fileListA = dir(obj.seqA);
                fileListA = fileListA(3:length(fileListA));
                fileListB = dir(obj.seqB);
                fileListB = fileListB(3:length(fileListB));
                fileNameA = fileListA(frameIndex).name;
                fileNameB = fileListB(frameIndex).name;
                IA = imread([obj.seqA fileNameA]);
                IB = imread([obj.seqB fileNameB]);
                
                valid = squeeze(obj.CP(frameIndex, obj.validCP(frameIndex, :) == 1, :));
                IA = DrawFeature(IA, valid(:, 1:2), [255 255 255], 'o');
                IB = DrawFeature(IB, valid(:, 3:4), [0 255 255], 's');
                warpA = obj.render1('a', IA, frameIndex, obj.Pa, obj.Ca, 1);
                warpB = obj.render1('b', IB, frameIndex, obj.Pb, obj.Cb, 1);
                
                imwrite(uint8(warpA), [outPath '/A' int2str(frameIndex) '.png']);
                imwrite(uint8(warpB), [outPath '/B' int2str(frameIndex) '.png']);
                warpA = uint8(warpA); warpB = uint8(warpB);
                warp = obj.imageBlending(warpA, warpB);
                maska = warpA(:, :, 1) == 0 & warpA(:, :, 2) == 0 & warpA(:, :, 3) == 0;
                maskb = warpB(:, :, 1) == 0 & warpB(:, :, 2) == 0 & warpB(:, :, 3) == 0;
                mask = maska & maskb;
                mask = repmat(mask, [1 1 3]);
                one = 255 * ones(size(mask));
                warp(mask) = one(mask);
                
                imwrite(warp, [outPath '/' int2str(frameIndex) '.png']);
            end
        end
        
        function CP_ = getStitchedCP(obj, head, tail)
            CP_ = zeros(tail - head + 1, obj.maxppf, 4);
            if obj.gap == 0
                obj.gap = 500;
            end
            for frameIndex = head:tail
                for k = 1:obj.maxppf
                    if obj.validCP(frameIndex, k) == 1
                        pa = obj.CP(frameIndex, k, 1:2);
                        pb = obj.CP(frameIndex, k, 3:4);
                        rowa = ceil(pa(2) / obj.quadHeight);
                        cola = ceil(pa(1) / obj.quadWidth);
                        rowb = ceil(pb(2) / obj.quadHeight);
                        colb = ceil(pb(1) / obj.quadWidth);
                        [paX, paY] = obj.transform(pa, obj.Offset * squeeze(obj.Pa(frameIndex, rowa, cola, :, :)) / squeeze(obj.Ca(frameIndex, rowa, cola, :, :)));
                        [pbX, pbY] = obj.transform(pb, obj.Offset * obj.H * squeeze(obj.Pb(frameIndex, rowb, colb, :, :)) / squeeze(obj.Cb(frameIndex, rowb, colb, :, :)));                        
                        CP_(frameIndex - head + 1, k, 1:2) = [paX + obj.gap paY + obj.gap];
                        CP_(frameIndex - head + 1, k, 3:4) = [pbX + obj.gap pbY + obj.gap];
                    else
                        CP_(frameIndex - head + 1, k, :) = -1; 
                    end
                end
            end
        end
        
        function addgrid(obj, outPathA, outPathB)
            if ~exist(outPathA, 'dir')
                mkdir(outPathA);
            end
            if ~exist(outPathB, 'dir')
                mkdir(outPathB);
            end
            width = 2;
            for frameIndex = 1:obj.nFrames %parfor
                disp(['rendering: # ' int2str(frameIndex)]);

                fileListA = dir(obj.seqA);
                fileListA = fileListA(3:length(fileListA));
                fileListB = dir(obj.seqB);
                fileListB = fileListB(3:length(fileListB));
                fileNameA = fileListA(frameIndex).name;
                fileNameB = fileListB(frameIndex).name;
                IA = imread([obj.seqA fileNameA]);
                IB = imread([obj.seqB fileNameB]);

                colorA = [133 206 244];
                colorB = [235 235 122];
                for row = 1:obj.meshSize
                    for col = 1:obj.meshSize
                        IA = obj.draw1Grid(IA, colorA, row, col, width); 
                        IB = obj.draw1Grid(IB, colorB, row, col, width);
                    end
                end
                imwrite(IA, [outPathA '/' sprintf('%03d',frameIndex) '.png']);
                imwrite(IB, [outPathB '/' sprintf('%03d',frameIndex) '.png']);
            end
        end
        
        function Igrid = draw1Grid(obj, I, color, row, col, width)
            Igrid = I;
            Igrid((row-1) * obj.quadHeight + 1 : (row-1) * obj.quadHeight + width, (col-1) * obj.quadWidth+1: col*obj.quadWidth, 1) = color(1);
            Igrid((row-1) * obj.quadHeight + 1 : (row-1) * obj.quadHeight + width, (col-1) * obj.quadWidth+1: col*obj.quadWidth, 2) = color(2);
            Igrid((row-1) * obj.quadHeight + 1 : (row-1) * obj.quadHeight + width, (col-1) * obj.quadWidth+1: col*obj.quadWidth, 3) = color(3);
            Igrid(row * obj.quadHeight - width + 1 : row * obj.quadHeight, (col-1) * obj.quadWidth + 1: col*obj.quadWidth, 1) = color(1);
            Igrid(row * obj.quadHeight - width + 1 : row * obj.quadHeight, (col-1) * obj.quadWidth + 1: col*obj.quadWidth, 2) = color(2);
            Igrid(row * obj.quadHeight - width + 1 : row * obj.quadHeight, (col-1) * obj.quadWidth + 1: col*obj.quadWidth, 3) = color(3);
            Igrid((row-1) * obj.quadHeight+1:row * obj.quadHeight, (col-1) * obj.quadWidth + 1:(col-1) * obj.quadWidth + width, 1) = color(1);
            Igrid((row-1) * obj.quadHeight+1:row * obj.quadHeight, (col-1) * obj.quadWidth + 1:(col-1) * obj.quadWidth + width, 2) = color(2);
            Igrid((row-1) * obj.quadHeight+1:row * obj.quadHeight, (col-1) * obj.quadWidth + 1:(col-1) * obj.quadWidth + width, 3) = color(3);
            Igrid((row-1) * obj.quadHeight+1:row * obj.quadHeight, col * obj.quadWidth - width + 1:col * obj.quadWidth, 1) = color(1);
            Igrid((row-1) * obj.quadHeight+1:row * obj.quadHeight, col * obj.quadWidth - width + 1:col * obj.quadWidth, 2) = color(2);
            Igrid((row-1) * obj.quadHeight+1:row * obj.quadHeight, col * obj.quadWidth - width + 1:col * obj.quadWidth, 3) = color(3);
            
            if row == 1
                Igrid((row-1) * obj.quadHeight + 1 : (row-1) * obj.quadHeight + width * 2, (col-1) * obj.quadWidth+1: col*obj.quadWidth, 1) = color(1);
                Igrid((row-1) * obj.quadHeight + 1 : (row-1) * obj.quadHeight + width * 2, (col-1) * obj.quadWidth+1: col*obj.quadWidth, 2) = color(2);
                Igrid((row-1) * obj.quadHeight + 1 : (row-1) * obj.quadHeight + width * 2, (col-1) * obj.quadWidth+1: col*obj.quadWidth, 3) = color(3);
            end
            
            if row == obj.meshSize
                Igrid(row * obj.quadHeight - 2 * width + 1 : row * obj.quadHeight, (col-1) * obj.quadWidth + 1: col*obj.quadWidth, 1) = color(1);
                Igrid(row * obj.quadHeight - 2 * width + 1 : row * obj.quadHeight, (col-1) * obj.quadWidth + 1: col*obj.quadWidth, 2) = color(2);
                Igrid(row * obj.quadHeight - 2 * width + 1 : row * obj.quadHeight, (col-1) * obj.quadWidth + 1: col*obj.quadWidth, 3) = color(3);
            end
            
            if col == 1
                Igrid((row-1) * obj.quadHeight+1:row * obj.quadHeight, (col-1) * obj.quadWidth + 1:(col-1) * obj.quadWidth + 2 * width, 1) = color(1);
                Igrid((row-1) * obj.quadHeight+1:row * obj.quadHeight, (col-1) * obj.quadWidth + 1:(col-1) * obj.quadWidth + 2 * width, 2) = color(2);
                Igrid((row-1) * obj.quadHeight+1:row * obj.quadHeight, (col-1) * obj.quadWidth + 1:(col-1) * obj.quadWidth + 2 * width, 3) = color(3);
            end
            
            if col == obj.meshSize
                Igrid((row-1) * obj.quadHeight+1:row * obj.quadHeight, col * obj.quadWidth - 2 * width + 1:col * obj.quadWidth, 1) = color(1);
                Igrid((row-1) * obj.quadHeight+1:row * obj.quadHeight, col * obj.quadWidth - 2 * width + 1:col * obj.quadWidth, 2) = color(2);
                Igrid((row-1) * obj.quadHeight+1:row * obj.quadHeight, col * obj.quadWidth - 2 * width + 1:col * obj.quadWidth, 3) = color(3);
            end
        end
        
        function imwarp = render1(obj, ab, I, frameIndex, P, C, drawGrid)
            src = Mesh(obj.videoHeight, obj.videoWidth, obj.quadWidth, obj.quadHeight);
            des = Mesh(obj.videoHeight, obj.videoWidth, obj.quadWidth, obj.quadHeight);
            %
            if ab == 'a'
                HH = obj.Offset;
            else
                HH = obj.Offset * obj.H;
            end
            
            % padding offset
            % ==============
            padding = 0;
            padOffset = zeros(obj.meshSize, obj.meshSize, 3, 3);
            for row = 1:obj.meshSize
                for col = 1:obj.meshSize
%                     padOffset(row, col, :, :) = squeeze(obj.Pa(1, row, col, :, :)) / squeeze(obj.Pa(padding, row, col, :, :));
                    padOffset(row, col, :, :) = eye(3);
                    if frameIndex > padding && frameIndex <= padding * 2
                        padOffset(row, col, :, :) = squeeze(padOffset(row, col, :, :)) * (squeeze(obj.Pa(ceil((frameIndex - padding)/1), row, col, :, :)) / squeeze(obj.Pa(padding, row, col, :, :)));
                    end
                    if frameIndex > obj.nFrames - padding * 2 && frameIndex <= obj.nFrames - padding
                        padOffset(row, col, :, :) = squeeze(padOffset(row, col, :, :)) * (squeeze(obj.Pb(obj.nFrames - ceil((obj.nFrames - padding - frameIndex)/1), row, col, :, :)) / squeeze(obj.Pb(obj.nFrames - padding, row, col, :, :)));
                    end
                end
            end
            avgOffset = zeros(obj.meshSize, obj.meshSize, 3, 3);
            for row = 1:obj.meshSize
                for col = 1:obj.meshSize
                    avgOffset(row, col, :, :) = avgOffset(row, col, :, :) + padOffset(row, col, :, :); 
                end
            end
            for row = 1:obj.meshSize
                for col = 1:obj.meshSize
                    avgOffset(row, col, :, :) = squeeze(avgOffset(row, col, :, :)) ./(obj.meshSize * obj.meshSize); 
                end
            end
            
            % ==============
            for i = 0 : obj.meshSize
                for j = 0 : obj.meshSize
                    x = i * obj.quadHeight + 1;
                    y = j * obj.quadWidth + 1;
                    if i == 0 && j == 0
                        B11 = squeeze(avgOffset(i+1, j+1, :, :)) * HH * squeeze(P(frameIndex, i+1, j+1, :, :)) ...
            / squeeze(C(frameIndex, i+1, j+1, :, :));
                        [xx11, yy11] = obj.transform([y x], B11);
                        des.setVertex(i,j, myPoint(xx11, yy11));
                        continue;
                    end
                    if i == 0 && j == obj.meshSize
                        B10 = squeeze(avgOffset(i+1, j, :, :)) * HH * squeeze(P(frameIndex, i+1, j, :, :))...
            / squeeze(C(frameIndex, i+1, j, :, :)) ;
                        %B10 = B10;

                        [xx10, yy10] = obj.transform([y x], B10);
                        des.setVertex(i,j, myPoint(xx10, yy10));
                        continue;
                    end
                    if i == 0 && j > 0 && j < obj.meshSize
                        B11 = squeeze(avgOffset(i+1, j+1, :, :)) * HH * squeeze(P(frameIndex, i+1, j+1, :, :))  ...
            / squeeze(C(frameIndex, i+1, j+1, :, :));
                        B10 = squeeze(avgOffset(i+1, j, :, :)) * HH * squeeze(P(frameIndex, i+1, j, :, :))  ...
            / squeeze(C(frameIndex, i+1, j, :, :));
                        [xx11, yy11] = obj.transform([y x], B11);
                        [xx10, yy10] = obj.transform([y x], B10);
                        des.setVertex(i,j, myPoint(mean([xx10 xx11]), mean([yy10 yy11])));
                        continue;
                    end
                    if i>0 && i < obj.meshSize && j == 0
                        B11 = squeeze(avgOffset(i+1, j+1, :, :)) * HH * squeeze(P(frameIndex, i+1, j+1, :, :))  ...
            / squeeze(C(frameIndex, i+1, j+1, :, :));
                        B01 = squeeze(avgOffset(i, j+1, :, :)) * HH * squeeze(P(frameIndex, i, j+1, :, :))  ...
            / squeeze(C(frameIndex, i, j+1, :, :));
                        [xx11, yy11] = obj.transform([y x], B11);
                        [xx01, yy01] = obj.transform([y x], B01);
                        des.setVertex(i,j, myPoint(mean([xx01 xx11]), mean([yy01 yy11])));
                        continue;
                    end
                    if i > 0 && i < obj.meshSize && j > 0 && j < obj.meshSize
                        B11 = squeeze(avgOffset(i+1, j+1, :, :)) * HH * squeeze(P(frameIndex, i+1, j+1, :, :))  ...
            / squeeze(C(frameIndex, i+1, j+1, :, :));
                        B01 = squeeze(avgOffset(i, j+1, :, :)) * HH * squeeze(P(frameIndex, i, j+1, :, :))  ...
            / squeeze(C(frameIndex, i, j+1, :, :));
                        B00 = squeeze(avgOffset(i, j, :, :)) * HH * squeeze(P(frameIndex, i, j, :, :))  ...
            / squeeze(C(frameIndex, i, j, :, :));
                        B10 = squeeze(avgOffset(i+1, j, :, :)) * HH * squeeze(P(frameIndex, i+1, j, :, :))  ...
            / squeeze(C(frameIndex, i+1, j, :, :));
                        [xx11, yy11] = obj.transform([y x], B11);
                        [xx01, yy01] = obj.transform([y x], B01);
                        [xx10, yy10] = obj.transform([y x], B10);
                        [xx00, yy00] = obj.transform([y x], B00);
                        [xx, yy] = obj.mergepoints([xx01, xx11, xx10, xx00], [yy01, yy11, yy10, yy00]);
                        des.setVertex(i,j, myPoint(xx, yy));
                        continue;
                    end
                    if i>0 && i < obj.meshSize && j == obj.meshSize
                        B00 = squeeze(avgOffset(i, j, :, :)) * HH * squeeze(P(frameIndex, i, j, :, :))  ...
            / squeeze(C(frameIndex, i, j, :, :));
                        B10 = squeeze(avgOffset(i+1, j, :, :)) * HH * squeeze(P(frameIndex, i+1, j, :, :))  ...
            / squeeze(C(frameIndex, i+1, j, :, :));
                        [xx10, yy10] = obj.transform([y x], B10);
                        [xx00, yy00] = obj.transform([y x], B00);
                        des.setVertex(i,j, myPoint(mean([xx10 xx00]), mean([yy10 yy00])));
                        continue;
                    end
                    if i == obj.meshSize && j == 0
                        B01 = squeeze(avgOffset(i, j+1, :, :)) * HH * squeeze(P(frameIndex, i, j+1, :, :))  ...
            / squeeze(C(frameIndex, i, j+1, :, :));
                        [xx01, yy01] = obj.transform([y x], B01);
                        des.setVertex(i,j, myPoint(xx01, yy01));
                        continue;
                    end
                    if i == obj.meshSize && j>0 && j < obj.meshSize
                        B00 = squeeze(avgOffset(i, j, :, :)) *  HH * squeeze(P(frameIndex, i, j, :, :))  ...
            / squeeze(C(frameIndex, i, j, :, :));
                        B01 = squeeze(avgOffset(i, j+1, :, :)) * HH * squeeze(P(frameIndex, i, j+1, :, :))  ...
            / squeeze(C(frameIndex, i, j+1, :, :));
                        [xx00, yy00] = obj.transform([y x], B00);
                        [xx01, yy01] = obj.transform([y x], B01);
                        des.setVertex(i,j, myPoint(mean([xx01 xx00]), mean([yy01 yy00])));
                        continue;
                    end
                    if i == obj.meshSize && j == obj.meshSize
                        B00 = squeeze(avgOffset(i, j, :, :)) * HH * squeeze(P(frameIndex, i, j, :, :)) ...
            / squeeze(C(frameIndex, i, j, :, :));
                        [xx00, yy00] = obj.transform([y x], B00); 
                        des.setVertex(i,j, myPoint(xx00, yy00));
                        continue;
                    end
                end
            end
            imwarp = zeros(obj.videoHeight+obj.gap*2,obj.videoWidth+obj.gap*2,3);

            for i=1:obj.meshSize
                for j=1:obj.meshSize
%                         disp([i,j]);
                    p0 = src.getVertex(i-1,j-1);
                    p1 = src.getVertex(i-1,j);
                    p2 = src.getVertex(i,j-1);
                    p3 = src.getVertex(i,j);

                    q0 = des.getVertex(i-1,j-1);
                    q1 = des.getVertex(i-1,j);
                    q2 = des.getVertex(i,j-1);
                    q3 = des.getVertex(i,j);

                    qd1 = Quad(p0,p1,p2,p3);
                    qd2 = Quad(q0,q1,q2,q3);
                    imwarp = quadWarp(obj,I,qd1,qd2, imwarp);
                end
            end 
            % grid
            if drawGrid
                imwarp = des.drawMesh(imwarp,obj.gap, ab);            
            end
        end
        
        function [xx, yy] = mergepoints(~, x, y)
            d = ones(4, 4) * 1e8;
            for i = 1:3
                for j = i+1:4
                    d(i, j) = (x(i) - x(j))^2 + (y(i) - y(j))^2;
                end
            end
            [minrow, min_i] = min(d);
            [~, min_j] = min(minrow);
            xx = 0.5 * (x(min_i(min_j))+ x(min_j));
            yy = 0.5 * (y(min_i(min_j))+ y(min_j));            
        end
        
        function imwarp = quadWarp(obj,im,q1,q2, imwarp)
            
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
            
            HH = homography_4pts(source',target');
            HH = HH./HH(3,3);

            imwarp = myWarp(minx,maxx,miny,maxy,double(im),imwarp,HH,obj.gap);
            imwarp = uint8(imwarp);
            
        end
        
        function output_canvas = imageBlending(~, warped_img1,warped_img2)
    
            w1 = imfill(im2bw(uint8(warped_img1), 0),'holes');
            w2 = imfill(im2bw(uint8(warped_img2), 0),'holes');

            w1 = mat2gray(w1);
            w2 = mat2gray(w2);

            warped_img1 = double(warped_img1);
            warped_img2 = double(warped_img2);
            output_canvas(:,:,1) = ((warped_img1(:,:,1).*w1)+(warped_img2(:,:,1).*w2))./(w1+w2);
            output_canvas(:,:,2) = ((warped_img1(:,:,2).*w1)+(warped_img2(:,:,2).*w2))./(w1+w2);
            output_canvas(:,:,3) = ((warped_img1(:,:,3).*w1)+(warped_img2(:,:,3).*w2))./(w1+w2);
            output_canvas = uint8(output_canvas);
            
        end
        
        function imwarp = warp1(obj, aorb, frameIndex, location_i, location_j, source, imwarp)
            if aorb == 'a'
                P = obj.Pa; 
                C = obj.Ca;
            else
                P = obj.Pb;
                C = obj.Cb;
            end
            B = squeeze(P(frameIndex, location_i, location_j, :, :)) / squeeze(C(frameIndex, location_i, location_j, :, :));
            if aorb == 'b'
                B = obj.H * B; 
            end
            B = obj.Offset * B;
            %B = B ./B(3,3);
            minx = (location_i - 1)*obj.quadHeight+1;
            maxx = minx + obj.quadHeight - 1;
            miny = (location_j - 1)*obj.quadWidth+1;
            maxy = miny + obj.quadWidth - 1;
            [x00, y00 ] = obj.transform([miny, minx], B);
            [x01, y01 ] = obj.transform([miny, maxx], B);
            [x10, y10 ] = obj.transform([maxy, minx], B);
            [x11, y11 ] = obj.transform([maxy, maxx], B);
            minx = min(x00, x01);minx = min(minx, x10);minx = min(minx, x11);
            miny = min(y00, y10);miny = min(miny, y01);miny = min(miny, y11);
            maxx = max(x10, x11);maxx = max(maxx, x01);maxx = max(maxx, x00);
            maxy = max(y01, y11);maxy = max(maxy, y10);maxy = max(maxy, y00);
            % imwarp = myWarp(minx, maxx, miny, maxy, double(source), obj.imwarp(frameIndex).cdata, inv(B), obj.gap);
            imwarp = myWarp(minx, maxx, miny, maxy, double(source), imwarp, inv(B), obj.gap);
        end
   
        
        
        function calcOmega(obj)
%             disp('computing omega...')
            for i = 1:obj.meshSize
                for j = 1:obj.meshSize
                    for t = 1:obj.nFrames                        
                        for r = t-obj.span:t+obj.span
%                             if t > obj.span && t <= obj.nFrames - obj.span
                            if r > 0 && r <= obj.nFrames
                                dPa = abs(obj.Ca(t,i,j,1,3) - obj.Ca(r,i,j,1,3)) + abs(obj.Ca(t,i,j,2,3) - obj.Ca(r,i,j,2,3));
                                dPb = abs(obj.Cb(t,i,j,1,3) - obj.Cb(r,i,j,1,3)) + abs(obj.Cb(t,i,j,2,3) - obj.Cb(r,i,j,2,3));
                                obj.w_a(t,r,i,j) = gaussmf(abs(t-r), [10 0]) * gaussmf(dPa, [400 0]);
                                obj.w_a(t,t,i,j) = 0;
                                obj.w_b(t,r,i,j) = gaussmf(abs(t-r), [10 0]) * gaussmf(dPb, [400 0]);
                                obj.w_b(t,t,i,j) = 0;
                            end
                        end
                        obj.gamma_a(t,i,j) = sum(obj.w_a(t,:,i,j)) * 2 * obj.smoothness;
                        obj.gamma_b(t,i,j) = sum(obj.w_b(t,:,i,j)) * 2 * obj.smoothness;
                        obj.gamma_a(t,i,j) = obj.gamma_a(t,i,j) + 1;
                        obj.gamma_b(t,i,j) = obj.gamma_b(t,i,j) + 1;
                        if (1 == obj.meshSize) || obj.stableW == 0
                            continue;
                        end
                        if ((i==1)||(i==obj.meshSize))&&((j==1)||(j==obj.meshSize))
                            obj.gamma_a(t,i,j) = obj.gamma_a(t,i,j) + 2 * 3 * obj.stableW;
                            obj.gamma_b(t,i,j) = obj.gamma_b(t,i,j) + 2 * 3 * obj.stableW;
                        end
                        if ((i==1)||(i==obj.meshSize))&&((j>1)&&(j<obj.meshSize))
                            obj.gamma_a(t,i,j) = obj.gamma_a(t,i,j) + 2 * 5 * obj.stableW;
                            obj.gamma_b(t,i,j) = obj.gamma_b(t,i,j) + 2 * 5 * obj.stableW;
                        end
                        if ((i>1)&&(i<obj.meshSize))&&((j==1)||(j==obj.meshSize))
                            obj.gamma_a(t,i,j) = obj.gamma_a(t,i,j) + 2 * 5 * obj.stableW;
                            obj.gamma_b(t,i,j) = obj.gamma_b(t,i,j) + 2 * 5 * obj.stableW;
                        end
                        if ((i>1)&&(i<obj.meshSize))&&((j>1)&&(j<obj.meshSize))
                            obj.gamma_a(t,i,j) = obj.gamma_a(t,i,j) + 2 * 8 * obj.stableW;
                            obj.gamma_b(t,i,j) = obj.gamma_b(t,i,j) + 2 * 8 * obj.stableW;
                        end
                    end
                end
            end            
        end
        
        function [x,y] = transform(~, xxyy, B)
            xx = xxyy(1); yy = xxyy(2);
            res = B * [xx;yy;1];
            x = res(1)/res(3);
            y = res(2)/res(3);
        end
        
        function CaPa = timesCa(obj, Pa)
            %
            x = Pa(1); y = Pa(2); frameIndex = Pa(3);
            if (x * y) ~= 0
                i = floor((x - 0.0001) / obj.quadWidth) + 1;
                j = floor((y - 0.0001) / obj.quadHeight) + 1;
                [x, y] = obj.transform(Pa, squeeze(obj.Ca_inv(frameIndex, i, j, :, :)));
            end
            CaPa = ones(3, 1);
            CaPa(1) = x; CaPa(2) = y;
        end
        
        function CaPa = timesCb(obj, Pa)
            %
            x = Pa(1); y = Pa(2); frameIndex = Pa(3);
            if (x * y) ~= 0
                i = floor((x - 0.0001) / obj.quadWidth) + 1;
                j = floor((y - 0.0001) / obj.quadHeight + 1);
                [x, y] = obj.transform(Pa, squeeze(obj.Cb_inv(frameIndex, i, j, :, :)));
            end
            CaPa = ones(3, 1);
            CaPa(1) = x; CaPa(2) = y;
        end
        
        
        
    end
end

