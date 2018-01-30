
classdef Quad < handle
    %QUAD Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    V00,V01,V10,V11;
    end
    
    methods
        function obj =Quad(V00,V01,V10,V11)
            obj.V00 = V00;
            obj.V01 = V01;
            obj.V10 = V10;
            obj.V11 = V11;
        end
        function a = isPointIn(obj,pt)
           in1 = isPointInTriangular(pt,obj.V00,obj.V01,obj.V11);
           in2 = isPointInTriangular(pt,obj.V00,obj.V10,obj.V11);
           if in1==1 || in2==1
               a = 1;
           else
               a = 0;    
           end
        end
        
        function a = isPointsIn(obj,pts)
            in1 = isPointsInTriangular(pts(:,1),pts(:,2),obj.V00,obj.V01,obj.V11);
            in2 = isPointsInTriangular(pts(:,1),pts(:,2),obj.V00,obj.V10,obj.V11);
            
            a = double(in1|in2);
        
        end
        
        function [coefficients,bool] = getBilinearCoordinates(obj,pt)
            	

            a_x = obj.V00.x - obj.V01.x - obj.V10.x + obj.V11.x;
            b_x = -obj.V00.x + obj.V01.x;
            c_x = -obj.V00.x + obj.V10.x;
            d_x = obj.V00.x - pt.x;
            
            a_y = obj.V00.y - obj.V01.y - obj.V10.y + obj.V11.y;
            b_y = -obj.V00.y + obj.V01.y;
            c_y = -obj.V00.y + obj.V10.y;
            d_y = obj.V00.y - pt.y;
            
            bigA = -a_y*b_x + b_y*a_x;
            bigB = -a_y*d_x - c_y*b_x + d_y*a_x +b_y*c_x;
            bigC = -c_y*d_x + d_y*c_x;
            
            tmp1 = -1;
            tmp2 = -1;
            tmp3 = -1;
            tmp4 = -1;

            if bigB*bigB - 4*bigA*bigC >= 0.0
                if abs(bigA) >= 0.000001
                    tmp1 = ( -bigB + sqrt(bigB*bigB - 4*bigA*bigC) ) / ( 2*bigA );
                    tmp2 = ( -bigB - sqrt(bigB*bigB - 4*bigA*bigC) ) / ( 2*bigA );
                else
                    tmp1 = -bigC/bigB;
                end
                
                if ( tmp1 >= -0.999999 && tmp1 <= 1.000001)
                    
                    tmp3 = -(b_y*tmp1 + d_y) / (a_y*tmp1 + c_y);
                    tmp4 = -(b_x*tmp1 + d_x) / (a_x*tmp1 + c_x);
                    if tmp3 >= -0.999999 && tmp3 <= 1.000001
                        k1 = tmp1;
                        k2 = tmp3;
                    else if tmp4 >= -0.999999 && tmp4 <= 1.000001
                            k1 = tmp1;
                            k2 = tmp4;
                        end
                    end
                end
                if ( tmp2 >= -0.999999 && tmp2 <= 1.000001)
                    
                    if tmp3 >= -0.999999 && tmp3 <= 1.000001
                        k1 = tmp2;
                        k2 = tmp3;
                    else if tmp4 >= -0.999999 && tmp4 <= 1.000001
                            k1 = tmp2;
                            k2 = tmp4;
                        end
                    end
                end
            end
	
            if k1>=-0.999999 && k1<=1.000001 && k2>=-0.999999 && k2<=1.000001
                
                coe1 = (1.0-k1)*(1.0-k2);
                coe2 = k1*(1.0-k2);
                coe3 = (1.0-k1)*k2;
                coe4 = k1*k2;
                
                coefficients(1) = coe1;
                coefficients(2) = coe2;
                coefficients(3) = coe3;
                coefficients(4) = coe4;
                
                bool = 1;
            else
                bool = 0;
            end
        
        end
        
        function minx = getMinX(obj)
            minx = min(obj.V00.x,obj.V01.x);
            minx = min(minx, obj.V10.x);
            minx = min(minx,obj.V11.x);
        end
        function maxx = getMaxX(obj)
            maxx = max(obj.V00.x,obj.V01.x);
            maxx = max(maxx,obj.V10.x);
            maxx = max(maxx,obj.V11.x);
        end
        function miny = getMinY(obj)
            miny = min(obj.V00.y,obj.V01.y);
            miny = min(miny,obj.V10.y);
            miny = min(miny,obj.V11.y);
        end
        function maxy = getMaxY(obj)
           maxy = max(obj.V00.y,obj.V01.y);
           maxy = max(maxy,obj.V10.y);
           maxy = max(maxy,obj.V11.y);
        end
        
        
    end
end


function a = isPointInTriangular(pt,V0,V1,V2)
    lambda1 = ((V1.y-V2.y)*(pt.x-V2.x) + (V2.x-V1.x)*(pt.y-V2.y)) / ((V1.y-V2.y)*(V0.x-V2.x) + (V2.x-V1.x)*(V0.y-V2.y));
	lambda2 = ((V2.y-V0.y)*(pt.x-V2.x) + (V0.x-V2.x)*(pt.y-V2.y)) / ((V2.y-V0.y)*(V1.x-V2.x) + (V0.x-V2.x)*(V1.y-V2.y));
    lambda3 = 1-lambda1-lambda2;
    if lambda1 >= 0.0 && lambda1 <= 1.0 && lambda2 >= 0.0 && lambda2 <= 1.0 && lambda3 >= 0.0 && lambda3 <= 1.0
        a = 1; %true
    else
        a = 0; %false
    end
end

function a = isPointsInTriangular(ptx,pty,V0,V1,V2)
    y12 = V1.y - V2.y;
    x21 = V2.x - V1.x;
    x02 = V0.x - V2.x;
    x21 = V2.x - V1.x;
    y02 = V0.y - V2.y;
    y20 = V2.y - V0.y;
    x12 = V1.x - V2.x;
    s1 = y12*x02 + x21*y02;
    s2 = y20*x12 + x02*y12;
    
    lambda1 = (y12.*(ptx-V2.x) + x21.*(pty-V2.y)) ./ s1;
    lambda2 = (y20.*(ptx-V2.x) + x02.*(pty-V2.y)) ./ s2;
    lambda3 = 1-lambda1-lambda2;
    
    condition1 = lambda1>=0.0;
    condition2 = lambda1<=1.0;
    condition3 = lambda2>=0.0;
    condition4 = lambda2<=1.0;
    condition5 = lambda3>=0.0;
    condition6 = lambda3<=1.0;
    
    a = condition1 & condition2 & condition3 & condition4 & condition5 & condition6;

end




