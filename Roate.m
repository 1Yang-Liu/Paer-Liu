function R = Roate(angleX,angleY,angleZ,S)
%%此为刚矩阵旋转矩阵,输入单位为°,矩阵为6×6
angleX=deg2rad(angleX);
angleY=deg2rad(angleY);
angleZ=deg2rad(angleZ);
Z=zeros(3);A=eye(3);
if S=='XYZ'
    Rx=[1,0,0;
        0,cos(angleX),-sin(angleX);
        0,sin(angleX),cos(angleX)];
    Ry=[cos(angleY),0,sin(angleY);
        0,1,0;
        -sin(angleY),0,cos(angleY)];
    Rz=[cos(angleZ),-sin(angleZ),0;
        sin(angleZ),cos(angleZ),0;
        0,0,1];
    R=[Rz*Ry*Rx,Z;Z,A];
else
    R=0;
end
end
