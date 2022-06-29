clc
clear
modu0 =[36.9792   36.9792   22.8933;
    7.0412    6.6603    6.6603;
    0.0973    0.3278    0.3278];

 moducf=[38.4118   38.4118   23.1805;
    7.1175    6.7067    6.7067;
    0.0952    0.3148    0.3148];
  
  modusf=[37.3039   37.3039   23.2252;
    7.1144    6.7338    6.7338;
    0.1011    0.3323    0.3323];
    
  moduts =[36.5998   36.5998   22.9225;
    7.0454    6.6762    6.6762;
    0.1027    0.3325    0.3325];
  modupn=[36.6881   36.6881   23.2996;
    7.0640    6.6988    6.6988;
    0.1024    0.3275    0.3275];

detts=moduts-modu0;
detcf=moducf-modu0;
detsf=modusf-modu0;
detpn=modupn-modu0;
DETE1=abs(detts(1,1))+abs(detcf(1,1))+abs(detsf(1,1))+abs(detpn(1,1));
DETE3=abs(detts(1,3))+abs(detcf(1,3))+abs(detsf(1,3))+abs(detpn(1,3));
DETG1=abs(detts(2,1))+abs(detcf(2,1))+abs(detsf(2,1))+abs(detpn(2,1));
DETG3=abs(detts(2,3))+abs(detcf(2,3))+abs(detsf(2,3))+abs(detpn(2,3));
E11=100*[detts(1,1),detcf(1,1),detsf(1,1),detpn(1,1)]/DETE1
E33=100*[detts(1,3),detcf(1,3),detsf(1,3),detpn(1,3)]/DETE3
G12=100*[detts(2,1),detcf(2,1),detsf(2,1),detpn(2,1)]/DETG1
G23=100*[detts(2,3),detcf(2,3),detsf(2,3),detpn(2,3)]/DETG3