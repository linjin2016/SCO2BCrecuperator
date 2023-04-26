clear
clc
%%������Ϊ��֤��ȷ����ѹ��SOC2���׶�ѭ��ϵͳ�ĳ��򣬶Ա�����ΪDyreby�Ĳ�ʿ����
%�������PCHE�������ṹΪֱͨ���������ɼ�֮ǰ��ֱͨ������ 
%Multi-objective optimization of supercritical carbon dioxode recompression
%Brayton cycle considering printed circuit recuperator design

%��֪���� ���ݵ���ѹ�������T2
intMC=0.89;%��ѹ��������Ч��
intRC=0.89;%��ѹ��������Ч��
intTB=0.93;%���ͻ��ĵ���Ч��
P1=7615;%��ѹ�������ѹ��7.69MPa
P2=25000;%��ѹ��������ѹ��20MPa
T1=32+273;%��ѹ��������¶�32��,K
S1=refpropm('S','T',T1,'P',P1,'co2');%��ѹ��������� J/kg K
Han1=refpropm('H','T',T1,'P',P1,'co2');%��ѹ��������� J/kg
Han2isen=refpropm('H','P',P2,'S',S1,'co2');%��ѹ�����س����ʣ�����ѹ������ K
Wcomp=(Han2isen-Han1)/intMC;
Han2=Han1+Wcomp;
T2=refpropm('T','P',P2,'H',Han2,'co2');
%����������������֪��Han�Ļ�������֪��ѹ��P���У�֪���¶�T��Han�޷����ѹ��P��������������
T6=700+273;%���ͻ�������¶ȣ�K
%QPHX=10000000;
FRC=0.3682;
%FRC=0.24;
ksaiLTR=0.95;
ksaiHTR=0.95;
GLTR=330; %���»���������������
GHTR=330; %���»���������������

%�����ֵ
P5=P2-500;
P9=P1+77;
P8=P9+30;
P3=P2;P4=P3;P10=P3;
P6=P5-245;
%T5=644;
P7=P1+200;

%  T8��ʼֵ���� ����������������
T8=470;
%  T8��ʼֵ���� ����������������

%���»�����
NumLTR=401;
twLTR=0.0005;
conwLTR=16.3;%��������ȵ���W/mK
DiamLTR=2.75;
Jh=12;Jc=12;
DUh=36;DUc=36; %��б�Ƕ�

ALTR=(pi*(DiamLTR/2)^2)/2;%ͨ�����������mm2
DLTR=pi*DiamLTR/2/(1+pi/2);%ͨ��ˮ��ֱ��,mm
AzhouchangLTR=pi*DiamLTR/2+DiamLTR;%ͨ���ܳ�,mm

mhLTR=GLTR*ALTR/1000000;%����ͨ���Ȳ�CO2����������
mcLTR=mhLTR*(1-FRC);%����ͨ�����CO2����������

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
LTRaveh=dPCLTR;%����0615

THLTR(NumLTR+1)=T8;
TCLTR(1)=T2;
PHLTR(NumLTR+1)=P8;
PHLTR(1)=P9;
PCLTR(1)=P2;
PCLTR(NumLTR+1)=P3;

%���»�����
NumHTR=401;
twHTR=0.0005;
conwHTR=16.3;%��������ȵ���W/mK
DiamHTR=2.75;
AHTR=(pi*(DiamHTR/2)^2)/2;%ͨ���������
DHTR=pi*DiamHTR/2/(1+pi/2);
AzhouchangHTR=pi*DiamHTR/2+DiamHTR;

mhHTR=GHTR*AHTR/1000000;%����ͨ���Ȳ�CO2����������
mcHTR=mhHTR;%����ͨ�����CO2����������

dPHHTR=zeros(1,NumHTR);
dPCHTR=zeros(1,NumHTR);
HanhHTR=zeros(1,NumHTR);
HancHTR=zeros(1,NumHTR);
LengthHTR=dPCHTR;
PHHTR=zeros(1,NumHTR+1);
PCHTR=zeros(1,NumHTR+1);
THHTR=PHHTR;TCHTR=PHHTR;
HTRaveh=dPCHTR;%����0615

%��ֵ����
    
%������ֵ����
mLTR=1;
nLTR=1;
pLTR=1;
mT8star=1;
nP5star=1;
nP7star=1;

    while abs((nP5star-P5)/P5)>0.00005|abs((nP7star-P7)/P7)>0.00005
     
    
     
     S6=refpropm('S','T',T6,'P',P6,'co2');%���ͻ������ J/kg K
     Han6=refpropm('H','T',T6,'P',P6,'co2');
     Han7isen=refpropm('H','P',P7,'S',S6,'co2');
     Wturb=intTB*(Han6-Han7isen);
     Han7=Han6-Wturb;
     T7=refpropm('T','P',P7,'H',Han7,'co2');%��ѹ�������¶ȣ�����ѹ������ K
     
    
     
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
    %�Ȳ�����
    CphLTR=refpropm('C','T',THLTR(i),'P',PHLTR(i),'co2');%�Ȳඨѹ������
    MiuhLTR=refpropm('V','T',THLTR(i),'P',PHLTR(i),'co2');%�Ȳද��ճ��
    DenhLTR=refpropm('D','T',THLTR(i),'P',PHLTR(i),'co2');%�Ȳ��ܶ�
    ConhLTR=refpropm('L','T',THLTR(i),'P',PHLTR(i),'co2');%�Ȳർ��ϵ��
    RehLTR=1000*mhLTR*DLTR/(MiuhLTR*ALTR);%�Ȳ���ŵ��
    ReynohL(i)=RehLTR;
    PrhLTR=refpropm('^','T',THLTR(i),'P',PHLTR(i),'co2');%�Ȳ���������
    HhLTR=refpropm('H','T',THLTR(i),'P',PHLTR(i),'co2');%�Ȳ�CO2��
    HanhLTR(i)=HhLTR;
    %������Բ���
    CpcLTR=refpropm('C','T',TCLTR(i),'P',PCLTR(i),'co2');%��ඨѹ������
    MiucLTR=refpropm('V','T',TCLTR(i),'P',PCLTR(i),'co2');%��ද��ճ��
    DencLTR=refpropm('D','T',TCLTR(i),'P',PCLTR(i),'co2');%����ܶ�
    ConcLTR=refpropm('L','T',TCLTR(i),'P',PCLTR(i),'co2');%��ർ��ϵ��
    RecLTR=1000*mcLTR*DLTR/(MiucLTR*ALTR);%�����ŵ��
    ReynocL(i)=RecLTR;
    PrcLTR=refpropm('^','T',TCLTR(i),'P',PCLTR(i),'co2');%�����������
    HcLTR=refpropm('H','T',TCLTR(i),'P',PCLTR(i),'co2');%���CO2�� 
    HancLTR(i)=HcLTR;
    
    %�Ȳ�Nu��f
    fhLTR=0.1924*(RehLTR^(-0.091)); 
    NuhLTR=0.184*(RehLTR^0.629)*(PrhLTR^0.317);
    %���Nu��f
    fcLTR=0.1924*(RecLTR^(-0.091)); 
    NucLTR=0.184*(RecLTR^0.629)*(PrcLTR^0.317);
    NucLTRjuzhen(i)=NucLTR;
    
    %����ϵ��
    hhLTR=1000*NuhLTR*ConhLTR/DLTR;%�Ȳഫ��ϵ��
    hcLTR=1000*NucLTR*ConcLTR/DLTR;%��ഫ��ϵ��
    hallLTR=1/(1/hhLTR+1/hcLTR+twLTR/conwLTR);%���廻��ϵ��
    LTRaveh(i)=hallLTR;%����0615
    
    delTQLTR=QLTR/NumLTR;%J
    delTLLTR=delTQLTR/hallLTR/(THLTR(i)-TCLTR(i))/(AzhouchangLTR/1000);%m
    LengthLTR(i)=delTLLTR;% m
    
    %�Ȳ�ѹ��
    deltaPhLTR=fhLTR*2*(LengthLTR(i))*mhLTR^2/(DenhLTR*(DLTR/1000)*(ALTR/1000000)^2)/1000; %����f��ʽ����ڵ��ѹ��,������1000��ת��Ϊkpa
    dPHLTR(i)=deltaPhLTR;%�˴�Ϊ�ڵ��ۼ�ѹ��
    PHLTR(i+1)=PHLTR(i)+deltaPhLTR;%������һ���ڵ��ѹ��
    %���ѹ��
    deltaPcLTR=fcLTR*2*LengthLTR(i)*i*mcLTR^2/(DencLTR*(DLTR/1000)*(ALTR/1000000)^2)/1000;%����f��ʽ����ڵ��ѹ��,,������1000��ת��Ϊkpa
    dPCLTR(i)=deltaPcLTR;%�˴�Ϊ���ڵ�ѹ��
    PCLTR(i+1)=PCLTR(1)-deltaPcLTR;%������һ���ڵ��ѹ��
    
    %%������ֵ����ڵ����(����ڵ��¶ȷ�����)
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
    
    %������ѹ��ģ�ͼ���T10
   P10=P3;
   S9=refpropm('S','T',T9,'P',P9,'co2');%���ͻ������ J/kg K
   Han9=refpropm('H','T',T9,'P',P9,'co2');%��ѹ��������� J/kg
   Han10isen=refpropm('H','P',P10,'S',S9,'co2');%��ѹ�����س����ʣ�����ѹ������ K
   Wrcomp=(Han10isen-Han9)/intRC;
   Han10=Han9+Wrcomp;
   T10=refpropm('T','P',P10,'H',Han10,'co2');
   
   %����T3��T10����T4
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
    
    %���»�����
    for i=1:NumHTR     
    %�Ȳ�����
    CphHTR=refpropm('C','T',THHTR(i),'P',PHHTR(i),'co2');%�Ȳඨѹ������
    MiuhHTR=refpropm('V','T',THHTR(i),'P',PHHTR(i),'co2');%�Ȳද��ճ��
    DenhHTR=refpropm('D','T',THHTR(i),'P',PHHTR(i),'co2');%�Ȳ��ܶ�
    ConhHTR=refpropm('L','T',THHTR(i),'P',PHHTR(i),'co2');%�Ȳർ��ϵ��
    RehHTR=1000*mhHTR*DHTR/(MiuhHTR*AHTR);%�Ȳ���ŵ��
    ReynohH(i)=RehHTR;
    PrhHTR=refpropm('^','T',THHTR(i),'P',PHHTR(i),'co2');%�Ȳ���������
    HhHTR=refpropm('H','T',THHTR(i),'P',PHHTR(i),'co2');%�Ȳ�CO2��
    HanhHTR(i)=HhHTR;
    %������Բ���
    CpcHTR=refpropm('C','T',TCHTR(i),'P',PCHTR(i),'co2');%��ඨѹ������
    MiucHTR=refpropm('V','T',TCHTR(i),'P',PCHTR(i),'co2');%��ද��ճ��
    DencHTR=refpropm('D','T',TCHTR(i),'P',PCHTR(i),'co2');%����ܶ�
    ConcHTR=refpropm('L','T',TCHTR(i),'P',PCHTR(i),'co2');%��ർ��ϵ��
    RecHTR=1000*mcHTR*DHTR/(MiucHTR*AHTR);%�����ŵ��
    ReynocH(i)=RecHTR;
    PrcHTR=refpropm('^','T',TCHTR(i),'P',PCHTR(i),'co2');%�����������
    HcHTR=refpropm('H','T',TCHTR(i),'P',PCHTR(i),'co2');%���CO2�� 
    HancHTR(i)=HcHTR;
    %�Ȳ�Nu��f
    fhHTR=0.1924*(RehHTR^(-0.091));%�Ȳ�f
    NuhHTR=0.184*(RehHTR^0.629)*(PrhHTR^0.317);%�Ȳ�Nu 
    %���Nu��f
    fcHTR=0.1924*(RecHTR^(-0.091)); %���f
    NucHTR=0.184*(RecHTR^0.629)*(PrcHTR^0.317);%���Nu
    NucHTRjuzhen(i)=NucHTR;
    
    %����ϵ��
    hhHTR=1000*NuhHTR*ConhHTR/DHTR;%�Ȳഫ��ϵ��
    hcHTR=1000*NucHTR*ConcHTR/DHTR;%��ഫ��ϵ��
    hallHTR=1/(1/hhHTR+1/hcHTR+twHTR/conwHTR);%���廻��ϵ��
    HTRaveh(i)=hallHTR;%����0615
    
    delTQHTR=QHTR/NumHTR;%J
    delTLHTR=delTQHTR/hallHTR/(THHTR(i)-TCHTR(i))/(AzhouchangHTR/1000);%m
    LengthHTR(i)=delTLHTR;% m
    
    %�Ȳ�ѹ��
    deltaPhHTR=fhHTR*2*LengthHTR(i)*mhHTR^2/(DenhHTR*(DHTR/1000)*(AHTR/1000000)^2)/1000; %����f��ʽ����ڵ��ѹ��
    dPHHTR(i)=deltaPhHTR;%�˴�Ϊ�ڵ��ۼ�ѹ��
    PHHTR(i+1)=PHHTR(i)+deltaPhHTR;%������һ���ڵ��ѹ��
    %���ѹ��
    deltaPcHTR=fcHTR*2*LengthHTR(i)*i*mcHTR^2/(DencHTR*DHTR/1000*(AHTR/1000000)^2)/1000;%����f��ʽ����ڵ��ѹ��
    dPCHTR(i)=deltaPcHTR;%�˴�Ϊ���ڵ�ѹ��
    PCHTR(i+1)=PCHTR(1)-deltaPcHTR;%������һ���ڵ��ѹ��
    
    %%������ֵ����ڵ����(����ڵ��¶ȷ�����)
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
   
   T0=293.15;%�����¶�  
   P0=101.33;%����ѹ�� 
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
   
   %�ɱ�����
   %���ͻ�ѹ�����ɱ�
   WTturb=mCO2*Wturb;%���ͻ��ܹ�
   WTcomp=mCO2*(1-FRC)*Wcomp+mCO2*FRC*Wrcomp;%ѹ�����ܹ�
   Costcomp=6898*(WTcomp/1000)^(0.7865);%ѹ�����ܳɱ�
   Costturb=7790*(WTturb/1000)^(0.6842);%���ͻ��ܳɱ�
   %��Դ��Դ�ɱ�
   UAPHE=QPHX/(T6-T5);%
   UAPC=QPC/(T9-T1);%
   CostPHE=3500*(UAPHE/1000);
   CostPC=2300*(UAPC/1000);
   %�������ɱ�
   %���»�����
   LengthL=sum(LengthLTR(:));
   Wch=3.3;
   Hch=4;
   ncnLTR=mCO2/mhLTR;
   VolumeLTR=Wch*Hch*LengthL*ncnLTR/1000000;
   CostLTR=7940*120*VolumeLTR;
   %���»�����
   LengthH=sum(LengthHTR(:));
   Wch=3.3;
   Hch=4;
   ncnHTR=mCO2/mhHTR;
   VolumeHTR=Wch*Hch*LengthH*ncnHTR/1000000;
   CostHTR=7940*120*VolumeHTR;
   CostRecup=CostLTR+CostHTR;
   %�ܳɱ�
   CostTotol=CostRecup+CostPC+CostPHE+Costcomp+Costturb;
  
     %����Ч��
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
   
   %����ϵͳ����
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
   
   %Ein0=QPHX*(1-T0/Treyuan);%��Դ�������
   Wnet=WTturb-WTcomp;%�������
   %intEx=Wnet/Ein0;%����Ч��
   
   Ex1=((Han1-Han0)-T0*(Shang1-Shang0));%1����� %ע�����Ϊ��λ����������ֵ
   Ex2=((Han2-Han0)-T0*(Shang2-Shang0));%2�����
   Ex3=((Han3-Han0)-T0*(Shang3-Shang0));%3�����
   Ex4=((Han4-Han0)-T0*(Shang4-Shang0));%4�����
   Ex5=((Han5-Han0)-T0*(Shang5-Shang0));%5�����
   Ex6=((Han6-Han0)-T0*(Shang6-Shang0));%6�����
   Ex7=((Han7-Han0)-T0*(Shang7-Shang0));%7�����
   Ex8=((Han8-Han0)-T0*(Shang8-Shang0));%8�����
   Ex9=((Han9-Han0)-T0*(Shang9-Shang0));%9�����
   Ex10=((Han10-Han0)-T0*(Shang10-Shang0));%10�����
   
   
   
   Exsunturbine=mCO2*(Ex6-Ex7)-WTturb;%���ֻ�������
   AAA=mCO2*(Ex6-Ex7);
   ExsunHTR=mCO2*(Ex4+Ex7-Ex5-Ex8);%���»�����������
   ExsunLTR=mCO2*(1-FRC)*Ex2+mCO2*Ex8-mCO2*(1-FRC)*Ex3-mCO2*Ex9;%���»�����������
  % ExsunPC=-QPC+mCO2*(1-FRC)*(Ex9-Ex1);
   ExsunmainCom=mCO2*(1-FRC)*Wcomp-mCO2*(1-FRC)*(Ex2-Ex1);%��ѹ����������
   ExsunreCom=mCO2*FRC*Wrcomp-mCO2*FRC*(Ex10-Ex9);%��ѹ����������
   %ExsunPHE=QPHX+mCO2*(Ex5-Ex6);
  
   ExsunPHE=(1-1/0.8)*(Ex5-Ex6);%������Դ������
   ExsunPC=(1-0.5)*(Ex9-Ex1);%��ȴ��������
   Exsunbujian=Exsunturbine+ExsunHTR+ExsunLTR+ExsunmainCom+ExsunreCom+ExsunPC+ExsunPHE;
   shuchu=[Exsunbujian,Exsunturbine,ExsunHTR,ExsunLTR,ExsunmainCom,ExsunreCom,ExsunPC,ExsunPHE];
   
   scresult=[ExsunHTR/1000000,ExsunLTR/1000000,Exsunturbine/1000000,ExsunmainCom/1000000,ExsunreCom/1000000,ExsunPHE/1000000,ExsunPC/1000000];
   
   result=[ksaiHTR,CostRecup/1000000,Costcomp/1000000,Costturb/1000000,CostPHE/1000000,CostPC/1000000];
   
    %����0615
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