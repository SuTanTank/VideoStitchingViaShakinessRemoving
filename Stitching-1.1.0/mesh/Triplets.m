classdef Triplets < handle
    %TRIPLET Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = public)
    nmax,m,n,nnz;
    nzval,ri,ci;
    end
    
    methods (Access = public)
        function obj=Triplets()
            obj.nmax = 0; 
            obj.m = 0;
            obj.n = 0;
            obj.nnz = 0;
            obj.ri = 0;
            obj.ci = 0;
        end
        function clear(obj)
            obj.nmax = 0;
            obj.m = 0;
            obj.n = 0;
            obj.nnz = 0;
            obj.ri = 0;
            obj.ci = 0;
        end
        function reserve(obj,rows,cols,num)
           clear(obj);
           obj.m = rows;
           obj.n = cols;
           obj.nmax = num;
           obj.nzval = zeros(num);
           obj.ri = zeros(num);
           obj.ci = zeros(num);
        end
        
        function add_row(obj,r,c,val)
           obj.ri(obj.nnz+1) = r;
           obj.ci(obj.nnz+1) = c;
           obj.nzval(obj.nnz+1) = val;
           obj.nnz = obj.nnz + 1;
           obj.nnz
           
        end
    end
end

