clc
clear
%定义材料参数
    Ef1=200; Ef3=30; Gf12=9; Gf31=11.881; vf12=0.20;vf23=0.27; %纤维参数
    Em=15.9; Gm=5.7; vm=0.395; %基体参数
    Ep=2.5*10^-5; Gp=1*10^-5; vp=0.25; %孔隙参数
    Dn=1.45; lf=800;df=7;    %针刺区域半径   %纤维的长细
%RVE的结构参数  Ln=9.07; pn=100/(Ln^2);  %论文Xu
tu=0.94; ts=0.42;  pn=1.2156*2; Ln=(100/pn)^0.5;  %论文Xu
Vcf=0.23;Vcp=0.09; Vcm=1-Vcf-Vcp; %论文XU
Vsf=0.10;Vsp=0.2;  Vsm=1-Vsf-Vsp;  %论文XU
Vnf=0.75;Vnp=0;    Vnm=1-Vnf-Vnp;  %论文XU

% Vnf=0.01;Vnp=0;Vnm=1-Vnf-Vnp;  %论文Meng
% tu=0.42;ts=0.42;pn=15.38;Ln=(100/pn)^0.5;   %论文Meng
% Vcf=0.2531;Vcp=0.09;Vcm=1-Vcf-Vcp; %论文Meng
% Vsf=0.1097;Vsp=0.2;Vsm=1-Vsf-Vsp;  %论文Meng

 %%%%%Chamis经验公式用于单向布模量计算

            Eem=(Em*Vcm+Ep*Vcp)/(Vcp+Vcm);vem=vm;Gem=Eem/[2*(1+vem)];  %等效基体计算
            Eu11=Vcf*Ef1+(1-Vcf)*Eem; Eu33=Eem/[1-Vcf^0.5*(1-Eem/Ef3)];Eu22=Eu33;
            Gu12=Gem/[1-Vcf^0.5*(1-Gem/Gf12)]; Gu23=Gem/[1-Vcf^0.5*(1-Gem/Gf31)];Gu13=Gu12;
            vu12=Vcf*vf12+(1-Vcf)*vem; vu23=Eu33/(2*Gu23)-1;vu13=vu12;
             Su=[1/Eu11,-vu12/Eu11,-vu13/Eu11,0,0,0;-vu12/Eu11,1/Eu22,-vu23/Eu22,0,0,0;
                 -vu13/Eu11,-vu23/Eu22,1/Eu33,0,0,0;
                 0,0,0,1/Gu23,0,0;0,0,0,0,1/Gu13,0;0,0,0,0,0,1/Gu12];%定义单向布的柔度矩阵
            UDcloth=[Eu11,Eu22,Eu33;Gu12,Gu23,Gu13;vu12,vu23,vu13];
      
%%%%Halpin-Tsai模型用于短切纤维毡模量计算

            Eem=(Em*Vsm+Ep*Vsp)/(Vsp+Vsm);vem=vm;Gem=Eem/[2*(1+vem)];
            n1=[Ef1/Eem-1]/[Ef1/Eem+2*(lf/df)];
            n2=(Ef3/Eem-1)/(Ef3/Eem+2);
            E1=Eem*[1+(2*lf/df)*n1*Vsf]/(1-n1*Vsf);
            E2=Eem*[1+2*n2*Vsf]/(1-n2*Vsf);
            Es=3*E1/8+5*E2/8;Gs=E1/8+E2/4;vs=Es/(2*Gs)-1;
            Ss=[1/Es,-vs/Es,-vs/Es,0,0,0;-vs/Es,1/Es,-vs/Es,0,0,0;-vs/Es,-vs/Es,1/Es,0,0,0;
                0,0,0,1/Gs,0,0;0,0,0,0,1/Gs,0;0,0,0,0,0,1/Gs];%定义纤维毡的柔度矩阵
            felt=[Es,Gs,vs];

%Chamis经验公式用于针刺区域模量计算

            Eem=(Em*Vnm+Ep*Vnp)/(Vnp+Vnm);vem=vm;Gem=Eem/[2*(1+vem)];  %等效基体计算
            En11=Vnf*Ef1+(1-Vnf)*Eem; En33=Eem/[1-Vnf^0.5*(1-Eem/Ef3)];En22=En33;
            Gn12=Gem/[1-Vnf^0.5*(1-Gem/Gf12)]; Gn23=Gem/[1-Vnf^0.5*(1-Gem/Gf31)];Gn13=Gn12;
            vn12=Vnf*vf12+(1-Vnf)*vem; vn23=En33/(2*Gn23)-1;vn13=vn12;
            Sn=[1/En33,-vn23/En22,-vn13/En11,0,0,0;-vn23/En22,1/En22,-vn12/En11,0,0,0;
                 -vn13/En11,-vn12/En11,1/En11,0,0,0;
                 0,0,0,1/Gn12,0,0;0,0,0,0,1/Gn13,0;0,0,0,0,0,1/Gn23];
            needle=[En11,En22,En33;Gn12,Gn23,Gn13;vn12,vn23,vn13];
           
%%%%经典层合板理论用于RVE模型模量计算
        An=3.1415*Dn^2/4;AR=Ln^2;   %针刺区域面积

%     %%[0/w/45/w/90/w/-45/w/0]n层合板角度  %论文Meng
%         Cu90=Roate(0,0,90,'XYZ')'*inv(Su)*Roate(0,0,90,'XYZ');%转角公式 %论文Meng
%         Cu45=Roate(0,0,45,'XYZ')'*inv(Su)*Roate(0,0,45,'XYZ');
%         Cum45=Roate(0,0,-45,'XYZ')'*inv(Su)*Roate(0,0,-45,'XYZ');
%         CU=[inv(Su)+Cu90+Cu45+Cum45]/4;  SU=inv(CU);    %转角后的单向布刚度矩阵。
   
    %%%[0/w/90/w]n层合板角度  %论文XU

          Cu90=Roate(0,0,90,'XYZ')'*inv(Su)*Roate(0,0,90,'XYZ');%转角公式 %论文Xu1
           CU=[inv(Su)+Cu90]/2;  SU=inv(CU);               %转角后的单向布刚度矩阵。
% %        SU=Su;CU=inv(SU);                                          %论文Xu2,单向，不需要转角
        CL=[tu/(tu+ts)]*CU+[ts/(tu+ts)]*inv(Ss);        %Laminite：刚度矩阵叠加。
        SR=(1-An/AR)*inv(CL)+[An/AR]*Sn;                %RVE：柔度矩阵叠加。
        ER11=1/SR(1,1); ER22=1/SR(2,2); ER33=1/SR(3,3);      %提取柔度矩阵中的模量和泊松比。
        GR12=1/SR(4,4); GR23=1/SR(5,5); GR31=1/SR(6,6);      
        vR12=-SR(2,1)*ER11; vR23=-SR(3,2)*ER22;vR23=-SR(1,3)*ER33;
modulus=[ER11,ER22,ER33;GR12,GR23,GR31;vR12,vR23,vR23];

modulus1=[(ER11+ER22)/2,(ER11+ER22)/2,ER33;GR12,(GR23+GR31)/2,(GR23+GR31)/2;vR12,vR23,vR23]
        CR=inv(SR);     %提取刚度矩阵中的刚度系数。
stiffness=[CR(1,1),CR(2,2),CR(3,3)];
