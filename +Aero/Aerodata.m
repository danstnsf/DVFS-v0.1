classdef Aerodata < handle

    properties
        %static aerodynamic coefficients (function of alpha and beta)
        CD
        CS
        CL
        Cl
        Cm
        Cn
        %dynamic aerodynamic coefficients derivatives (function of rotation rate)
        CDq
        CLq
        Cmq

        CSp
        Clp
        Cnp

        CSr
        Clr
        Cnr
        
        %surface control derivatives (function of control surface
        %deflection)
        CDd
        CSd
        CLd
        Cld
        Cmd
        Cnd
    end
    
    methods
        function obj = Aerodata()
            if nargin
                
            end
        end
        function from_vsp_polar(obj,polarfile)
            if nargin    
            data=readtable(polarfile,"VariableNamingRule","preserve");
                alpha=unique(data.AoA)';
                beta=unique(data.Beta)';
                obj.CD=Aero.Coeff('CD',beta,alpha,[]);
                obj.CS=Aero.Coeff('CS',beta,alpha,[]);
                obj.CL=Aero.Coeff('CL',beta,alpha,[]);
                obj.Cl=Aero.Coeff('Cl',beta,alpha,[]);
                obj.Cm=Aero.Coeff('Cm',beta,alpha,[]);
                obj.Cn=Aero.Coeff('Cn',beta,alpha,[]);
                for i=1:length(beta)
                    obj.CL.bkvar1=beta;
                    obj.CL.bkvar2=alpha;
                    obj.CL.data(:,i)=groupfilter(data,'AoA',@(x) x==beta(i),"Beta").CL;
                    
                    obj.CD.bkvar1=beta;
                    obj.CD.bkvar2=alpha;
                    obj.CD.data(:,i)=groupfilter(data,'AoA',@(x) x==beta(i),"Beta").CDtot;
                    
                    obj.CS.bkvar1=beta;
                    obj.CS.bkvar2=alpha;
                    obj.CS.data(:,i)=groupfilter(data,'AoA',@(x) x==beta(i),"Beta").CS;
                    
                    obj.Cl.bkvar1=beta;
                    obj.Cl.bkvar2=alpha;
                    obj.Cl.data(:,i)=groupfilter(data,'AoA',@(x) x==beta(i),"Beta").CMl;
                    
                    obj.Cm.bkvar1=beta;
                    obj.Cm.bkvar2=alpha;
                    obj.Cm.data(:,i)=groupfilter(data,'AoA',@(x) x==beta(i),"Beta").CMm;
                    
                    obj.Cn.bkvar1=beta;
                    obj.Cn.bkvar2=alpha;
                    obj.Cn.data(:,i)=groupfilter(data,'AoA',@(x) x==beta(i),"Beta").CMn;
                end
                obj.CD.ft;
                obj.CS.ft;
                obj.CL.ft;
                obj.Cl.ft;
                obj.Cm.ft;
                obj.Cn.ft;
            end
        end
        function from_vsp_stab(obj,stabfile)
            if nargin               
                fileid=fopen(stabfile,'r');
                if fileid == -1
                    error('Não foi possível abrir o arquivo.');
                end
                content=fread(fileid,'*char')';
                fclose(fileid);
                lines=strsplit(content,'\n');
                linecont=cell(length(lines),1);
                for i=1:length(lines)
                    words=strsplit(lines{i});
                    linecont{i}=words;
                end
                
                n=(length(linecont)-1)/58;
                
                d=zeros(n,2);
                c=cell(n,1);
                for i=1:n
                    d(i,1)=str2double(cell2mat(linecont{11+58*(i-1)}(2)));
                    d(i,2)=str2double(cell2mat(linecont{12+58*(i-1)}(2)));
                    for j=1:12
                        temp=cell2table(linecont{(37+j)+58*(i-1)}(1:10));
                        c{i}=[c{i};temp];
                    end
                    c{i}.Properties.VariableNames=["Coef","Base","dalp","dbeta","dp","dq","dr","dM","dU","dctrl"];
                    c{i}=[array2table([d(i,1)*ones(12,1) d(i,2)*ones(12,1)],'VariableNames',{'alpha','beta'}) c{i}];
                    if i==1
                        L=c{i};
                    else
                        L=[L;c{i}];
                    end
                end
                L=L(:,[3,1,2,4,5,6,7,8,9,10,11,12]);
                L.Base=str2double(L.Base);
                L.dalp=str2double(L.dalp);
                L.dbeta=str2double(L.dbeta);
                L.dp=str2double(L.dp);
                L.dq=str2double(L.dq);
                L.dr=str2double(L.dr);
                L.dM=str2double(L.dM);
                L.dU=str2double(L.dU);
                L.dctrl=str2double(L.dctrl);
                
                MCD=groupfilter(L,"alpha",@(x) x=="CD","Coef");
                MCS=groupfilter(L,"alpha",@(x) x=="CS","Coef");
                MCL=groupfilter(L,"alpha",@(x) x=="CL","Coef");
                MCl=groupfilter(L,"alpha",@(x) x=="CMl","Coef");
                MCm=groupfilter(L,"alpha",@(x) x=="CMm","Coef");
                MCn=groupfilter(L,"alpha",@(x) x=="CMn","Coef");
                
                dyn_beta=unique(L.beta)';
                dyn_alpha=unique(L.alpha)';
                
                for i=1:length(dyn_beta)
                    temp=groupfilter(MCD,"alpha",@(x) x==dyn_beta(i),"beta");
                    for k=1:3
                        CDx(:,i,k)=table2array(temp(:,6+k));
                    end
                    CDx(:,i,4)=table2array(temp(:,12));
                end
                
                for i=1:length(dyn_beta)
                    temp=groupfilter(MCS,"alpha",@(x) x==dyn_beta(i),"beta");
                    for k=1:3
                        CSx(:,i,k)=table2array(temp(:,6+k));
                    end
                    CSx(:,i,4)=table2array(temp(:,12));
                end
                
                for i=1:length(dyn_beta)
                    temp=groupfilter(MCL,"alpha",@(x) x==dyn_beta(i),"beta");
                    for k=1:3
                        CLx(:,i,k)=table2array(temp(:,6+k));
                    end
                    CLx(:,i,4)=table2array(temp(:,12));
                end
                
                for i=1:length(dyn_beta)
                    temp=groupfilter(MCl,"alpha",@(x) x==dyn_beta(i),"beta");
                    for k=1:3
                        Clx(:,i,k)=table2array(temp(:,6+k));
                    end
                    Clx(:,i,4)=table2array(temp(:,12));
                end
                
                for i=1:length(dyn_beta)
                    temp=groupfilter(MCm,"alpha",@(x) x==dyn_beta(i),"beta");
                    for k=1:3
                        Cmx(:,i,k)=table2array(temp(:,6+k));
                    end
                    Cmx(:,i,4)=table2array(temp(:,12));
                end
                
                for i=1:length(dyn_beta)
                    temp=groupfilter(MCn,"alpha",@(x) x==dyn_beta(i),"beta");
                    for k=1:3
                        Cnx(:,i,k)=table2array(temp(:,6+k));
                    end
                    Cnx(:,i,4)=table2array(temp(:,12));
                end

                obj.CDq=Aero.Coeff('CDq',dyn_alpha,dyn_beta,CDx(:,:,2));
                obj.CDq.ft;
%                 obj.CDq.bkvar1=dyn_alpha;
%                 obj.CDq.bkvar2=dyn_beta;
%                 obj.CDq.data=CDx(:,:,2);
               
                obj.CSp=Aero.Coeff('CSp',dyn_alpha,dyn_beta,CSx(:,:,1));
                obj.CSp.ft;
%                 obj.CSp.bkvar1=dyn_alpha;
%                 obj.CSp.bkvar2=dyn_beta;
%                 obj.CSp.data=CSx(:,:,1);
                obj.CSr=Aero.Coeff('CSr',dyn_alpha,dyn_beta,CSx(:,:,3));
                obj.CSr.ft;
%                 obj.CSr.bkvar1=dyn_alpha;
%                 obj.CSr.bkvar2=dyn_beta;
%                 obj.CSr.data=CSx(:,:,3);                  
                
                obj.CLq=Aero.Coeff('CLq',dyn_alpha,dyn_beta,CLx(:,:,2));
                obj.CLq.ft;
%                 obj.CLq.bkvar1=dyn_alpha;
%                 obj.CLq.bkvar2=dyn_beta;
%                 obj.CLq.data=CLx(:,:,2);
                
                obj.Clp=Aero.Coeff('Clp',dyn_alpha,dyn_beta,Clx(:,:,1));
                obj.Clp.ft;
%                 obj.Clp.bkvar1=dyn_alpha;
%                 obj.Clp.bkvar2=dyn_beta;
%                 obj.Clp.data=Clx(:,:,1);
                
                obj.Clr=Aero.Coeff('Clr',dyn_alpha,dyn_beta,Clx(:,:,3));
                obj.Clr.ft;
%                 obj.Clr.bkvar1=dyn_alpha;
%                 obj.Clr.bkvar2=dyn_beta;
%                 obj.Clr.data=Clx(:,:,3);
                
                obj.Cmq=Aero.Coeff('Cmq',dyn_alpha,dyn_beta,Cmx(:,:,2));
                obj.Cmq.ft;
%                 obj.Cmq.bkvar1=dyn_alpha;
%                 obj.Cmq.bkvar2=dyn_beta;
%                 obj.Cmq.data=Cmx(:,:,2);
                
                obj.Cnp=Aero.Coeff('Cnp',dyn_alpha,dyn_beta,Cnx(:,:,1));
                obj.Cnp.ft;
%                 obj.Cnp.bkvar1=dyn_alpha;
%                 obj.Cnp.bkvar2=dyn_beta;
%                 obj.Cnp.data=Cnx(:,:,1);
                
                obj.Cnr=Aero.Coeff('Cnr',dyn_alpha,dyn_beta,Cnx(:,:,3));
                obj.Cnr.ft;
%                 obj.Cnr.bkvar1=dyn_alpha;
%                 obj.Cnr.bkvar2=dyn_beta;
%                 obj.Cnr.data=Cnx(:,:,3);
                
                obj.CDd=Aero.Coeff('CDd',dyn_alpha,dyn_beta,CDx(:,:,4));
                obj.CDd.ft;
%                 obj.CDd.bkvar1=dyn_alpha;
%                 obj.CDd.bkvar2=dyn_beta;
%                 obj.CDd.data=CDx(:,:,4);
                
                obj.CSd=Aero.Coeff('CSd',dyn_alpha,dyn_beta,CSx(:,:,4));
                obj.CSd.ft;
%                 obj.CSd.bkvar1=dyn_alpha;
%                 obj.CSd.bkvar2=dyn_beta;
%                 obj.CSd.data=CSx(:,:,4);
                
                obj.CLd=Aero.Coeff('CLd',dyn_alpha,dyn_beta,CLx(:,:,4));
                obj.CLd.ft;
%                 obj.CLd.bkvar1=dyn_alpha;
%                 obj.CLd.bkvar2=dyn_beta;
%                 obj.CLd.data=CLx(:,:,4);
                
                obj.Cld=Aero.Coeff('Cld',dyn_alpha,dyn_beta,Clx(:,:,4));
                obj.Cld.ft;
%                 obj.Cld.bkvar1=dyn_alpha;
%                 obj.Cld.bkvar2=dyn_beta;
%                 obj.Cld.data=Clx(:,:,4);
                
                obj.Cmd=Aero.Coeff('Cmd',dyn_alpha,dyn_beta,Cmx(:,:,4));
                obj.Cmd.ft;
%                 obj.Cmd.bkvar1=dyn_alpha;
%                 obj.Cmd.bkvar2=dyn_beta;
%                 obj.Cmd.data=Cmx(:,:,4);
                
                obj.Cnd=Aero.Coeff('Cnd',dyn_alpha,dyn_beta,Cnx(:,:,4));
                obj.Cnd.ft;
%                 obj.Cnd.bkvar1=dyn_alpha;
%                 obj.Cnd.bkvar2=dyn_beta;
%                 obj.Cnd.data=Cnx(:,:,4);
            end
        end

        function S=aerostate(obj,alpha,beta,w,V,c,b,d)
            p=w(1);
            q=w(2);
            r=w(3);
            
            qf=0.5*q*c/V;
            pf=0.5*p*b/V;
            rf=0.5*r*b/V;

            S=[obj.CD.fn(beta,alpha)+obj.CDq.fn(beta,alpha)*qf+obj.CDd.fn(beta,alpha)*d
               obj.CS.fn(beta,alpha)+obj.CSp.fn(beta,alpha)*pf+obj.CSr.fn(beta,alpha)*rf+obj.CSd.fn(beta,alpha)*d
               obj.CL.fn(beta,alpha)+obj.CLq.fn(beta,alpha)*qf+obj.CLd.fn(beta,alpha)*d
               obj.Cl.fn(beta,alpha)+obj.Clp.fn(beta,alpha)*pf+obj.Clr.fn(beta,alpha)*rf+obj.Cld.fn(beta,alpha)*d
               obj.Cm.fn(beta,alpha)+obj.Cmq.fn(beta,alpha)*qf+obj.Cmd.fn(beta,alpha)*d
               obj.Cn.fn(beta,alpha)+obj.Cnp.fn(beta,alpha)*pf+obj.Cnr.fn(beta,alpha)*rf+obj.Cnd.fn(beta,alpha)*d];
        end
    end
end

