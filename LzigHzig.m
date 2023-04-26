clear
clc
%%本程序为验证正确的再压缩SOC2布雷顿循环系统的程序，对比文献为Dyreby的博士论文
%本程序的PCHE换热器结构为直通道，参数可见之前的直通道文献 
%Multi-objective optimization of supercritical carbon dioxode recompression
%Brayton cycle considering printed circuit recuperator design

%已知参数 根据等熵压缩，算的T2
intMC=0.89;%主压缩机等熵效率
intRC=0.89;%再压缩机等熵效率
intTB=0.93;%膨胀机的等熵效率
P1=7615;%主压缩机入口压力7.69MPa
P2=25000;%主压缩机出口压力20MPa
T1=32+273;%主压缩机入口温度32度,K
S1=refpropm('S','T',T1,'P',P1,'co2');%主压缩机入口熵 J/kg K
Han1=refpropm('H','T',T1,'P',P1,'co2');%主压缩机入口焓 J/kg
Han2isen=refpropm('H','P',P2,'S',S1,'co2');%主压缩等熵出口焓，等熵压缩条件 K
Wcomp=(Han2isen-Han1)/intMC;
Han2=Han1+Wcomp;
T2=refpropm('T','P',P2,'H',Han2,'co2');
%！！！！！！！！知道Han的话，必须知道压力P才行，知道温度T和Han无法获得压力P！！！！！！！
T6=700+273;%膨胀机的入口温度，K
%QPHX=10000000;
FRC=0.3682;
%FRC=0.24;
ksaiLTR=0.95;
ksaiHTR=0.95;
GLTR=330; %低温回热器的质量流率
GHTR=330; %高温回热器的质量流率

%定义初值
P5=P2-500;
P9=P1+77;
P8=P9+30;
P3=P2;P4=P3;P10=P3;
P6=P5-245;
%T5=644;
P7=P1+200;

%  T8初始值定义 ！！！！！！！！
T8=470;
%  T8初始值定义 ！！！！！！！！

%低温回热器
NumLTR=401;
twLTR=0.0005;
conwLTR=16.3;%固体壁面热导率W/mK
DiamLTR=2.75;
Jh=12;Jc=12;
DUh=36;DUc=36; %倾斜角度

ALTR=(pi*(DiamLTR/2)^2)/2;%通道截面面积，mm2
DLTR=pi*DiamLTR/2/(1+pi/2);%通道水力直径,mm
AzhouchangLTR=pi*DiamLTR/2+DiamLTR;%通道周长,mm

mhLTR=GLTR*ALTR/1000000;%单个通道热侧CO2的质量流量
mcLTR=mhLTR*(1-FRC);%单个通道冷侧CO2的质量流量

dPHLTR=zeros(1,NumLTR);
dPCLTR=zeros(1,NumLTR);
HanhLTR=zeros(1,NumLTR);
HancLTR=zeros(1,NumLTR);

ReynohL=HancLTR;
ReynocL=HancLTR;
ReynohH=HancLTR;
ReynocH=HancLTR;

LengthLTR=dPCLTR;
PHLTR=zeros(1,NumLTR+1);
PCLTR=zeros(1,NumLTR+1);
THLTR=PHLTR;TCLTR=PHLTR;
NucLTRjuzhen=TCLTR;
NucHTRjuzhen=TCLTR;
LTRaveh=dPCLTR;%新增0615

THLTR(NumLTR+1)=T8;
TCLTR(1)=T2;
PHLTR(NumLTR+1)=P8;
PHLTR(1)=P9;
PCLTR(1)=P2;
PCLTR(NumLTR+1)=P3;

%高温回热器
NumHTR=401;
twHTR=0.0005;
conwHTR=16.3;%固体壁面热导率W/mK
DiamHTR=2.75;
AHTR=(pi*(DiamHTR/2)^2)/2;%通道截面面积
DHTR=pi*DiamHTR/2/(1+pi/2);
AzhouchangHTR=pi*DiamHTR/2+DiamHTR;

mhHTR=GHTR*AHTR/1000000;%单个通道热侧CO2的质量流量
mcHTR=mhHTR;%单个通道冷侧CO2的质量流量

dPHHTR=zeros(1,NumHTR);
dPCHTR=zeros(1,NumHTR);
HanhHTR=zeros(1,NumHTR);
HancHTR=zeros(1,NumHTR);
LengthHTR=dPCHTR;
PHHTR=zeros(1,NumHTR+1);
PCHTR=zeros(1,NumHTR+1);
THHTR=PHHTR;TCHTR=PHHTR;
HTRaveh=dPCHTR;%新增0615

%初值定义
    
%收敛初值定义
mLTR=1;
nLTR=1;
pLTR=1;
mT8star=1;
nP5star=1;
nP7star=1;

    while abs((nP5star-P5)/P5)>0.00005|abs((nP7star-P7)/P7)>0.00005
     
    
     
     S6=refpropm('S','T',T6,'P',P6,'co2');%膨胀机入口熵 J/kg K
     Han6=refpropm('H','T',T6,'P',P6,'co2');
     Han7isen=refpropm('H','P',P7,'S',S6,'co2');
     Wturb=intTB*(Han6-Han7isen);
     Han7=Han6-Wturb;
     T7=refpropm('T','P',P7,'H',Han7,'co2');%主压缩出口温度，等熵压缩条件 K
     
    
     
    while abs((mT8star-T8)/T8)>0.00001
         
     THLTR(NumLTR+1)=T8;
     TCLTR(1)=T2;
     PHLTR(NumLTR+1)=P8;
     PHLTR(1)=P9;
     PCLTR(1)=P2;
     PCLTR(NumLTR+1)=P3;
        
     HanThinPHLTR=refpropm('H','T',THLTR(NumLTR+1),'P',PHLTR(NumLTR+1),'co2');
     HanTcinPHLTR=refpropm('H','T',TCLTR(1),'P',PHLTR(1),'co2');
     HanThinPCLTR=refpropm('H','T',THLTR(NumLTR+1),'P',PCLTR(NumLTR+1),'co2');
     HanTcinPCLTR=refpropm('H','T',TCLTR(1),'P',PCLTR(1),'co2');
     delHmaxLTR=min(mhLTR*(1-FRC)*(HanThinPCLTR-HanTcinPCLTR),mhLTR*(HanThinPHLTR-HanTcinPHLTR));
     QLTR=ksaiLTR*delHmaxLTR;     
     Han8=refpropm('H','T',T8,'P',P8,'co2');
     Han9=Han8-QLTR/mhLTR;
     T9=refpropm('T','P',P9,'H',Han9,'co2');
     THLTR(1)=T9;
   
    for i=1:NumLTR 
    %热侧物性
    CphLTR=refpropm('C','T',THLTR(i),'P',PHLTR(i),'co2');%热侧定压比热容
    MiuhLTR=refpropm('V','T',THLTR(i),'P',PHLTR(i),'co2');%热侧动力粘度
    DenhLTR=refpropm('D','T',THLTR(i),'P',PHLTR(i),'co2');%热侧密度
    ConhLTR=refpropm('L','T',THLTR(i),'P',PHLTR(i),'co2');%热侧导热系数
    RehLTR=1000*mhLTR*DLTR/(MiuhLTR*ALTR);%热侧雷诺数
    ReynohL(i)=RehLTR;
    PrhLTR=refpropm('^','T',THLTR(i),'P',PHLTR(i),'co2');%热侧普朗特数
    HhLTR=refpropm('H','T',THLTR(i),'P',PHLTR(i),'co2');%热侧CO2焓
    HanhLTR(i)=HhLTR;
    %冷侧物性参数
    CpcLTR=refpropm('C','T',TCLTR(i),'P',PCLTR(i),'co2');%冷侧定压比热容
    MiucLTR=refpropm('V','T',TCLTR(i),'P',PCLTR(i),'co2');%冷侧动力粘度
    DencLTR=refpropm('D','T',TCLTR(i),'P',PCLTR(i),'co2');%冷侧密度
    ConcLTR=refpropm('L','T',TCLTR(i),'P',PCLTR(i),'co2');%冷侧导热系数
    RecLTR=1000*mcLTR*DLTR/(MiucLTR*ALTR);%冷侧雷诺数
    ReynocL(i)=RecLTR;
    PrcLTR=refpropm('^','T',TCLTR(i),'P',PCLTR(i),'co2');%冷侧普朗特数
    HcLTR=refpropm('H','T',TCLTR(i),'P',PCLTR(i),'co2');%冷侧CO2焓 
    HancLTR(i)=HcLTR;
    
    %热侧Nu和f
    fhLTR=0.1924*(RehLTR^(-0.091)); 
    NuhLTR=0.184*(RehLTR^0.629)*(PrhLTR^0.317);
    %冷侧Nu和f
    fcLTR=0.1924*(RecLTR^(-0.091)); 
    NucLTR=0.184*(RecLTR^0.629)*(PrcLTR^0.317);
    NucLTRjuzhen(i)=NucLTR;
    
    %传热系数
    hhLTR=1000*NuhLTR*ConhLTR/DLTR;%热侧传热系数
    hcLTR=1000*NucLTR*ConcLTR/DLTR;%冷侧传热系数
    hallLTR=1/(1/hhLTR+1/hcLTR+twLTR/conwLTR);%整体换热系数
    LTRaveh(i)=hallLTR;%新增0615
    
    delTQLTR=QLTR/NumLTR;%J
    delTLLTR=delTQLTR/hallLTR/(THLTR(i)-TCLTR(i))/(AzhouchangLTR/1000);%m
    LengthLTR(i)=delTLLTR;% m
    
    %热侧压降
    deltaPhLTR=fhLTR*2*(LengthLTR(i))*mhLTR^2/(DenhLTR*(DLTR/1000)*(ALTR/1000000)^2)/1000; %根据f公式计算节点间压差,最后除以1000是转化为kpa
    dPHLTR(i)=deltaPhLTR;%此处为节点累计压差
    PHLTR(i+1)=PHLTR(i)+deltaPhLTR;%计算下一个节点的压力
    %冷侧压降
    deltaPcLTR=fcLTR*2*LengthLTR(i)*i*mcLTR^2/(DencLTR*(DLTR/1000)*(ALTR/1000000)^2)/1000;%根据f公式计算节点间压差,,最后除以1000是转化为kpa
    dPCLTR(i)=deltaPcLTR;%此处为各节点压差
    PCLTR(i+1)=PCLTR(1)-deltaPcLTR;%计算下一个节点的压力
    
    %%根据焓值计算节点参数(计算节点温度方法二)
    AFLTR=(AzhouchangLTR/1000)*LengthLTR(i);%m2
    HanhLTR(i+1)=HanhLTR(i)+hallLTR*AFLTR*(THLTR(i)-TCLTR(i))/mhLTR;
    HancLTR(i+1)=HancLTR(i)+hallLTR*AFLTR*(THLTR(i)-TCLTR(i))/mcLTR;
    THLTR(i+1)=THLTR(i)+(HanhLTR(i+1)-HanhLTR(i))/CphLTR;
    TCLTR(i+1)=TCLTR(i)+(HancLTR(i+1)-HancLTR(i))/CpcLTR;  
    end
    
    
    P8=PHLTR(NumLTR+1);
    %T8=THLTR(NumLTR+1);
    %T9=THLTR(1);
    T3=TCLTR(NumLTR+1);
    P3=PCLTR(NumLTR+1);
    
    %根据再压缩模型计算T10
   P10=P3;
   S9=refpropm('S','T',T9,'P',P9,'co2');%膨胀机入口熵 J/kg K
   Han9=refpropm('H','T',T9,'P',P9,'co2');%主压缩机入口焓 J/kg
   Han10isen=refpropm('H','P',P10,'S',S9,'co2');%主压缩等熵出口焓，等熵压缩条件 K
   Wrcomp=(Han10isen-Han9)/intRC;
   Han10=Han9+Wrcomp;
   T10=refpropm('T','P',P10,'H',Han10,'co2');
   
   %根据T3和T10计算T4
   P4=P3;
   Han3=refpropm('H','T',T3,'P',P3,'co2');
   Han4=(1-FRC)*Han3+FRC*Han10;
   T4=refpropm('T','P',P4,'H',Han4,'co2');
   
   THHTR(NumHTR+1)=T7;
   TCHTR(1)=T4;
   PHHTR(1)=P8;
   PCHTR(1)=P4;
   PCHTR(NumHTR+1)=P5;
   PHHTR(NumHTR+1)=P7;
   
   HanThinPHHTR=refpropm('H','T',THHTR(NumHTR+1),'P',PHHTR(NumHTR+1),'co2');
   HanTcinPHHTR=refpropm('H','T',TCHTR(1),'P',PHHTR(1),'co2');
   HanThinPCHTR=refpropm('H','T',THHTR(NumHTR+1),'P',PCHTR(NumHTR+1),'co2');
   HanTcinPCHTR=refpropm('H','T',TCHTR(1),'P',PCHTR(1),'co2');
   
   delHmaxHTR=min(mcHTR*(HanThinPCHTR-HanTcinPCHTR),mhHTR*(HanThinPHHTR-HanTcinPHHTR));
   QHTR=ksaiHTR*delHmaxHTR;     
   Han8=Han7-QHTR/mhHTR;
   T8star=refpropm('T','P',PHHTR(1),'H',Han8,'co2');
   
   mT8star=T8star;
   if abs((mT8star-T8)/T8)>0.00001
    aT8=T8star-T8;
    dqT8=abs(0.5*aT8);
    if aT8>0
        dT8=T8+dqT8;
    elseif aT8<0
        dT8=T8-dqT8;
    end
    T8=dT8;
   end
   %THLTR(1)=T8;
    end
    
    PHHTR(1)=P8;
    THHTR(1)=T8;
    PCHTR(1)=P4;
    TCHTR(1)=T4;
    
    %高温回热器
    for i=1:NumHTR     
    %热侧物性
    CphHTR=refpropm('C','T',THHTR(i),'P',PHHTR(i),'co2');%热侧定压比热容
    MiuhHTR=refpropm('V','T',THHTR(i),'P',PHHTR(i),'co2');%热侧动力粘度
    DenhHTR=refpropm('D','T',THHTR(i),'P',PHHTR(i),'co2');%热侧密度
    ConhHTR=refpropm('L','T',THHTR(i),'P',PHHTR(i),'co2');%热侧导热系数
    RehHTR=1000*mhHTR*DHTR/(MiuhHTR*AHTR);%热侧雷诺数
    ReynohH(i)=RehHTR;
    PrhHTR=refpropm('^','T',THHTR(i),'P',PHHTR(i),'co2');%热侧普朗特数
    HhHTR=refpropm('H','T',THHTR(i),'P',PHHTR(i),'co2');%热侧CO2焓
    HanhHTR(i)=HhHTR;
    %冷侧物性参数
    CpcHTR=refpropm('C','T',TCHTR(i),'P',PCHTR(i),'co2');%冷侧定压比热容
    MiucHTR=refpropm('V','T',TCHTR(i),'P',PCHTR(i),'co2');%冷侧动力粘度
    DencHTR=refpropm('D','T',TCHTR(i),'P',PCHTR(i),'co2');%冷侧密度
    ConcHTR=refpropm('L','T',TCHTR(i),'P',PCHTR(i),'co2');%冷侧导热系数
    RecHTR=1000*mcHTR*DHTR/(MiucHTR*AHTR);%冷侧雷诺数
    ReynocH(i)=RecHTR;
    PrcHTR=refpropm('^','T',TCHTR(i),'P',PCHTR(i),'co2');%冷侧普朗特数
    HcHTR=refpropm('H','T',TCHTR(i),'P',PCHTR(i),'co2');%冷侧CO2焓 
    HancHTR(i)=HcHTR;
    %热侧Nu和f
    fhHTR=0.1924*(RehHTR^(-0.091));%热侧f
    NuhHTR=0.184*(RehHTR^0.629)*(PrhHTR^0.317);%热侧Nu 
    %冷侧Nu和f
    fcHTR=0.1924*(RecHTR^(-0.091)); %冷侧f
    NucHTR=0.184*(RecHTR^0.629)*(PrcHTR^0.317);%冷侧Nu
    NucHTRjuzhen(i)=NucHTR;
    
    %传热系数
    hhHTR=1000*NuhHTR*ConhHTR/DHTR;%热侧传热系数
    hcHTR=1000*NucHTR*ConcHTR/DHTR;%冷侧传热系数
    hallHTR=1/(1/hhHTR+1/hcHTR+twHTR/conwHTR);%整体换热系数
    HTRaveh(i)=hallHTR;%新增0615
    
    delTQHTR=QHTR/NumHTR;%J
    delTLHTR=delTQHTR/hallHTR/(THHTR(i)-TCHTR(i))/(AzhouchangHTR/1000);%m
    LengthHTR(i)=delTLHTR;% m
    
    %热侧压降
    deltaPhHTR=fhHTR*2*LengthHTR(i)*mhHTR^2/(DenhHTR*(DHTR/1000)*(AHTR/1000000)^2)/1000; %根据f公式计算节点间压差
    dPHHTR(i)=deltaPhHTR;%此处为节点累计压差
    PHHTR(i+1)=PHHTR(i)+deltaPhHTR;%计算下一个节点的压力
    %冷侧压降
    deltaPcHTR=fcHTR*2*LengthHTR(i)*i*mcHTR^2/(DencHTR*DHTR/1000*(AHTR/1000000)^2)/1000;%根据f公式计算节点间压差
    dPCHTR(i)=deltaPcHTR;%此处为各节点压差
    PCHTR(i+1)=PCHTR(1)-deltaPcHTR;%计算下一个节点的压力
    
    %%根据焓值计算节点参数(计算节点温度方法二)
    AFHTR=(AzhouchangHTR/1000)*LengthHTR(i);
    HanhHTR(i+1)=HanhHTR(i)+hallHTR*AFHTR*(THHTR(i)-TCHTR(i))/mhHTR;
    HancHTR(i+1)=HancHTR(i)+hallHTR*AFHTR*(THHTR(i)-TCHTR(i))/mcHTR;
    THHTR(i+1)=THHTR(i)+(HanhHTR(i+1)-HanhHTR(i))/CphHTR;
    TCHTR(i+1)=TCHTR(i)+(HancHTR(i+1)-HancHTR(i))/CpcHTR;  
    end
    
    
    P7star=PHHTR(NumHTR+1);
    P5star=PCHTR(NumHTR+1);
    T5=TCHTR(NumHTR+1);
    T7=THHTR(NumHTR+1);
    
    nP5star=P5star;
    nP7star=P7star;
    
    if abs((nP5star-P5)/P5)>0.00005
    P5=nP5star;
    P6=P5-250;
    end
    
    if abs((nP7star-P7)/P7)>0.00005
    P7=nP7star;
    end
    
    end
   
   T0=293.15;%环境温度  
   P0=101.33;%环境压力 
   Han0=refpropm('H','T',T0,'P',P0,'co2');
   Han1=refpropm('H','T',T1,'P',P1,'co2');
   Han2=refpropm('H','T',T2,'P',P2,'co2');
   Han3=refpropm('H','T',T3,'P',P3,'co2');
   Han4=refpropm('H','T',T4,'P',P4,'co2');
    
   Han7=refpropm('H','T',T7,'P',P7,'co2');
   Han8=refpropm('H','T',T8,'P',P8,'co2');
   Han10=refpropm('H','T',T10,'P',P10,'co2');
    
   Han6=refpropm('H','T',T6,'P',P6,'co2');
   Han5=refpropm('H','T',T5,'P',P5,'co2');
   TMD=(Han6-Han5);
   mCO2=74.968;
   QPHX=mCO2*(Han6-Han5);
   %mCO2=QPHX/(Han6-Han5);
   
   Han9=refpropm('H','T',T9,'P',P9,'co2');
   QPC=mCO2*(1-FRC)*(Han9-Han1);
   
   intTH=1-(QPC/QPHX);
   
   %成本计算
   %膨胀机压缩机成本
   WTturb=mCO2*Wturb;%膨胀机总功
   WTcomp=mCO2*(1-FRC)*Wcomp+mCO2*FRC*Wrcomp;%压缩机总功
   Costcomp=6898*(WTcomp/1000)^(0.7865);%压缩机总成本
   Costturb=7790*(WTturb/1000)^(0.6842);%膨胀机总成本
   %热源冷源成本
   UAPHE=QPHX/(T6-T5);%
   UAPC=QPC/(T9-T1);%
   CostPHE=3500*(UAPHE/1000);
   CostPC=2300*(UAPC/1000);
   %回热器成本
   %低温回热器
   LengthL=sum(LengthLTR(:));
   Wch=3.3;
   Hch=4;
   ncnLTR=mCO2/mhLTR;
   VolumeLTR=Wch*Hch*LengthL*ncnLTR/1000000;
   CostLTR=7940*120*VolumeLTR;
   %高温回热器
   LengthH=sum(LengthHTR(:));
   Wch=3.3;
   Hch=4;
   ncnHTR=mCO2/mhHTR;
   VolumeHTR=Wch*Hch*LengthH*ncnHTR/1000000;
   CostHTR=7940*120*VolumeHTR;
   CostRecup=CostLTR+CostHTR;
   %总成本
   CostTotol=CostRecup+CostPC+CostPHE+Costcomp+Costturb;
  
     %火用效率
   Shang0=refpropm('S','T',T0,'P',P0,'co2');
   Shang1=refpropm('S','T',T1,'P',P1,'co2');
   Shang2=refpropm('S','T',T2,'P',P2,'co2');
   Shang3=refpropm('S','T',T3,'P',P3,'co2');
   Shang4=refpropm('S','T',T4,'P',P4,'co2');
    
   Shang7=refpropm('S','T',T7,'P',P7,'co2');
   Shang8=refpropm('S','T',T8,'P',P8,'co2');
   Shang9=refpropm('S','T',T9,'P',P9,'co2');
   Shang10=refpropm('S','T',T10,'P',P10,'co2');
    
   Shang6=refpropm('S','T',T6,'P',P6,'co2');
   Shang5=refpropm('S','T',T5,'P',P5,'co2');
   
   %计算系统熵增
%    szeng1=(1-FRC)*mCO2*(Shang2-Shang1);
%    szeng2=(1-FRC)*mCO2*(Shang3-Shang2);
%    szeng3=mCO2*Shang4-(1-FRC)*mCO2*Shang3-FRC*mCO2*Shang10;
%    szeng4=mCO2*(Shang5-Shang4);
%    szeng5=mCO2*(Shang6-Shang5);
%    szeng6=mCO2*(Shang7-Shang6);
%    szeng7=mCO2*(Shang8-Shang7);
%    szeng8=mCO2*(Shang9-Shang8);
%    szeng9=(1-FRC)*mCO2*(Shang1-Shang9);
%    szeng10=FRC*mCO2*(Shang10-Shang9);
%    szenghe=szeng1+szeng2+szeng3+szeng4+szeng5+szeng6++szeng7++szeng8++szeng9++szeng10;
   
   
  % Treyuan=1273;
   
   %Ein0=QPHX*(1-T0/Treyuan);%热源输入火用
   Wnet=WTturb-WTcomp;%净输出功
   %intEx=Wnet/Ein0;%火用效率
   
   Ex1=((Han1-Han0)-T0*(Shang1-Shang0));%1点火用 %注意此外为单位质量流量的值
   Ex2=((Han2-Han0)-T0*(Shang2-Shang0));%2点火用
   Ex3=((Han3-Han0)-T0*(Shang3-Shang0));%3点火用
   Ex4=((Han4-Han0)-T0*(Shang4-Shang0));%4点火用
   Ex5=((Han5-Han0)-T0*(Shang5-Shang0));%5点火用
   Ex6=((Han6-Han0)-T0*(Shang6-Shang0));%6点火用
   Ex7=((Han7-Han0)-T0*(Shang7-Shang0));%7点火用
   Ex8=((Han8-Han0)-T0*(Shang8-Shang0));%8点火用
   Ex9=((Han9-Han0)-T0*(Shang9-Shang0));%9点火用
   Ex10=((Han10-Han0)-T0*(Shang10-Shang0));%10点火用
   
   
   
   Exsunturbine=mCO2*(Ex6-Ex7)-WTturb;%汽轮机火用损
   AAA=mCO2*(Ex6-Ex7);
   ExsunHTR=mCO2*(Ex4+Ex7-Ex5-Ex8);%高温回热器火用损
   ExsunLTR=mCO2*(1-FRC)*Ex2+mCO2*Ex8-mCO2*(1-FRC)*Ex3-mCO2*Ex9;%低温回热器火用损
  % ExsunPC=-QPC+mCO2*(1-FRC)*(Ex9-Ex1);
   ExsunmainCom=mCO2*(1-FRC)*Wcomp-mCO2*(1-FRC)*(Ex2-Ex1);%主压缩机火用损
   ExsunreCom=mCO2*FRC*Wrcomp-mCO2*FRC*(Ex10-Ex9);%再压缩机火用损
   %ExsunPHE=QPHX+mCO2*(Ex5-Ex6);
  
   ExsunPHE=(1-1/0.8)*(Ex5-Ex6);%高温热源火用损
   ExsunPC=(1-0.5)*(Ex9-Ex1);%冷却器火用损
   Exsunbujian=Exsunturbine+ExsunHTR+ExsunLTR+ExsunmainCom+ExsunreCom+ExsunPC+ExsunPHE;
   shuchu=[Exsunbujian,Exsunturbine,ExsunHTR,ExsunLTR,ExsunmainCom,ExsunreCom,ExsunPC,ExsunPHE];
   
   scresult=[ExsunHTR/1000000,ExsunLTR/1000000,Exsunturbine/1000000,ExsunmainCom/1000000,ExsunreCom/1000000,ExsunPHE/1000000,ExsunPC/1000000];
   
   result=[ksaiHTR,CostRecup/1000000,Costcomp/1000000,Costturb/1000000,CostPHE/1000000,CostPC/1000000];
   
    %新增0615
   LTRcoldP=PCLTR(1)-PCLTR(NumLTR+1);
   LTRhotP=PHLTR(NumLTR+1)-PHLTR(1);
   HTRcoldP=PCHTR(1)-PCHTR(NumHTR+1);
   HTRhotP=PHHTR(NumHTR+1)-PHHTR(1);
   LTRaveK=mean(LTRaveh(:));
   HTRaveK=mean(HTRaveh(:));
   LTRtemcha=THLTR-TCLTR;
   LTRjiadian=min(LTRtemcha);
   HTRtemcha=THHTR-TCHTR;
   HTRjiadian=min(HTRtemcha);
   
   NucHTRavg=mean(NucHTRjuzhen);
   NucLTRavg=mean(NucLTRjuzhen);
   
   RESTS=[LTRcoldP,LTRhotP,HTRcoldP,HTRhotP,LTRaveK,HTRaveK,NucLTRavg,NucHTRavg,LTRjiadian,HTRjiadian,VolumeLTR,VolumeHTR,intTH*100,CostTotol/1000000,Exsunbujian/1000000];
   zongjieshuju=RESTS.';
   
   AME=[ksaiHTR,intTH*100,CostTotol/1000000,Exsunbujian/1000000];
   BME=[T1,T2,T3,T4,T5,T6,T7,T8,T9,T10];
   CME=[P1,P2,P3,P4,P5,P6,P7,P8,P9,P10];
   DME=[Ex1,Ex2,Ex3,Ex4,Ex5,Ex6,Ex7,Ex8,Ex9,Ex10];
   EME=[Shang1,Shang2,Shang3,Shang4,Shang5,Shang6,Shang7,Shang8,Shang9,Shang10];