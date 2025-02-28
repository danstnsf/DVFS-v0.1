classdef Coeff < handle
    % This is a aerodynamic coefficient class. It transforms your coefficient
    % data stored in 'datafile' (txt, csv, xls ...) to the combination of the
    % variable breakpoints (limited to two in this implementation) and a table
    % variable with the 'datafile' content.
    properties
        name
        bkvar1=[]
        bkvar2=[]
        data=[]
        
        ft
        fn
    end
    
    methods
        function obj = Coeff(name,breakpointsvar1,breakpointsvar2,data)
            if nargin
                obj.bkvar1 = breakpointsvar1;
                if nargin==4
                    obj.bkvar2 = breakpointsvar2;
                end
                obj.data=data;
                obj.name=name;
            end
        end

        function ft=get.ft(obj)
            if isempty(obj.bkvar2)
                ft=fit(obj.bkvar1,obj.data,'poly3');
            else
                [x,y]=meshgrid(obj.bkvar1,obj.bkvar2);
                ft=fit([x(:) y(:)],obj.data(:),'poly33');
                obj.fn=Aux.fitfcn(ft.p00,ft.p10,ft.p01,ft.p20,ft.p11,ft.p02,ft.p30,ft.p21,ft.p12,ft.p03);
            end
        end
    end
end

