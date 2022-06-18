clc;
clear;
addpath('Functions')
addpath('HalfCellData')
addpath('Input')

CellID=66824; %
HasCheckups=0; %Is 1 if the cell has checkup cycles and 0 if the cell has no checkup cycles
UpperorLower=1; %1 if fixed upper cutoff cells, 0 if fixed lower cutoff cells. 1 or 0 is fine for fixed middle cutoff cells. Keep the value at 1 if your cells are 100% DOD.
mode0=1; %0 for slippage/mass versus time using a 2D scan, 1 for Chi Square map, 2 for slippage/mass versus time using a 4D scan
ScanNumber=6; %6 %The number of total cycle numbers that are dV/dQ scanned. Each are approximately equally spaced.
Prec1=5; %Precision of active positive mass grid. Decrease number in order to decrease computational time.
Prec2=5; %Precision of positive slippage grid. Decrease number in order to decrease computational time.

maxV=4.1; %Maximum voltage in the data. Delete all possible noise above maxV.
minV=3.0; %Minimum voltage in the data. Delete all possible noise below minV.


%Load curves
A = readmatrix(sprintf('%d_QC.txt',CellID)); %Input capacity vs cycle curves 
B = readmatrix(sprintf('%d_VQ.txt',CellID)); %Input capacity vs voltage curves
C = readmatrix(sprintf('%d_IT.txt',CellID)); %Input current vs time curves

%{
%Load curves
A = readmatrix('66836_QC.txt'); %Input capacity vs cycle curves 
B = readmatrix('66836_VQ.txt'); %Input capacity vs voltage curves
C = readmatrix('66836_IT.txt'); %Input current vs time curves
%}

%Input negative curve
A0 = readmatrix('NG.csv');

%Input positive curve
B0 = readmatrix('622A.csv');

if mode0==0 || mode0==1
Negative_Mass=1.025; %Value that can be changed when in mode 0 or 1 (value is fixed) 
Negative_Slippage=-1; %Value that can be changed when in mode 0 or 1 (value is fixed) 

%---------------Ask Roby before changing anything below this line------------------
%Unless you know what you are doing.

O=Negative_Mass;
S=-Negative_Slippage;
SS=S;
OO=O;
else
end

if HasCheckups==1
%If the cell has checkup cycles   
D=Findcheckup(UpperorLower,A,B,C); %Function that find the checkup cycles. D=[n,T,Q,V,I];

else
%If the cell has no checkup cycles
Q=A(:,2); %Capacity
n=A(:,1); %Cycle number
V=B(:,2); %Voltage
I=[0;C(:,2)]; %Current
T=[0;C(:,1)]; %Time




D=[n,T,Q,V,I];
D=D(D(:,5)>0,:); %Only remove the discharge curves
end




%Clean data
D = D(all(D,2),:); %Remove rows that have zeros
D = D(D(:,4) <= maxV, :); %Remove all rows with V > maxV Volts
D = D(D(:,4) >= minV, :); %Remove all rows with V < minV Volts

[H,L]=size(D);

l=2;
i0=0;
%Identify each cycle number with a index l (letter l), each voltage within a cycle
%with an index k. Form a 3D matrix. The middle index of the 3D matrix correspond to a specific
%parameter 
for i=1:H-1

    diff=D(i+1,1)-D(i,1);  
    
 if diff>0.5
        l=l+1;
        i0=i;
    else
       k=i-i0; 
    end
    DD0(k,:,l)=D(i,:);
end 

LL=l;
DD0(DD0==0)=NaN; %Replace the zeros in the 3D matrix with an empty value

Sc=ScanNumber-1;
%Select only some cycles to scan instead of all of them. 
s=1;
    for l=2:floor(LL/Sc)-2:LL
      
    DD(:,:,s)=DD0(:,:,l);
    s=s+1;
end


%---------------------------------------------------------------------
clearvars -except DD mode0 A0 B0 Prec1 Prec2 CellID O S 

[H3,X,L3]=size(DD);

for s=1:L3

if mode0==0 || mode0==1 %|| s~=1    %add a not s==1 conditions (use Scan4D only for first cycle)...modify Scan2D accordingly
[MM,NN,GG,Pmass,Pslipp0]=Scan2D(DD,A0,B0,s,Prec1,Prec2,O,S);
OO=O;
SS=S;
else
[MM,NN,OO,SS,GG,Pmass,Pslipp0]=Scan4D(DD,A0,B0,s,Prec1,Prec2);  
S=SS;
O=OO;
end

[X,dV1,Xe,dVne] = dVdQSolution(DD,A0,B0,s,Prec1,Prec2,MM,NN,OO,SS);

  

PslippM(s)=-2.5*MM/Prec2; %Positive slippage
PmassN(s)=0.9+0.05*(NN-1)/Prec1; %Positive active mass
NmassO(s)=0.98+OO*0.02;
NslippS(s)=-SS;
Capacity(s)=max(DD(:,3,s));
cycle(s)=DD(1,1,s);
time(s)=DD(1,2,s)./3600;

Slipp(s)=-(PslippM(s)-NslippS(s));
%Slipp0=Slipp(1);
%Pmass0=PmassN(1);
%LiInv(s)=(Slipp(s)-Slipp0)+250.*(1-(PmassN(s)./Pmass0)); %Update 250 with a variable before sending to other students
%LiInv(s)=(Slipp(s)-28)+250.*(1-(PmassN(s)./1.39)); %Update 250 with a variable before sending to other students

if mode0==1
SlippageMap=-(Pslipp0-S); %Slippage for 2D map only
else
end    
cycle0=cycle(s);
time0=round(time(s));

fpath = 'Output';

figure(1)
 plot(X,dV1)
    xlim([8 250])
    hold on
   plot(Xe,dVne)
    xlim([8 250])
    ylim([0 0.025])
    xlabel('Capacity(mAh)')
    ylabel('dV/dQ')
    hold off 
  saveas(gcf,fullfile(fpath,sprintf('%d_dVdQ_4D_cycle = %d.jpg',CellID,cycle0)));

if mode0==1
LGG=-log(GG);

figure(2)
x0=500;
y0=200;
width=1080*1.5;
height=500*1.5;

set(gcf,'position',[x0,y0,width,height])
subplot(2,3,s);
image(Pmass,SlippageMap,LGG,'CDataMapping','scaled')

%xlim([1.4 1.5])
%ylim([-70 -40])
caxis([-0.5 5])
xlabel('Positive electrode mass (g)','FontSize',14) 
ylabel('Absolute slippage (mAh)','FontSize',14) 
set(gca,'FontSize',11)
set(gca,'YDir','normal')
c = colorbar;
c.Label.String = '-log(\chi^2)';
c.Label.FontSize=14;
colormap jet
title(sprintf('Cycle = %d; Time = %d h',cycle0,time0)) 
saveas(gcf,fullfile(fpath,sprintf('%d_matrix%02d_cycle=%d_time=%f.jpg',CellID,s,cycle0,time0)));
%pause(0.05)

else
end
end

Slipp0 = interp1(time,Slipp,0,'linear','extrap');
Pmass0 = interp1(time,PmassN,0,'linear','extrap');
%Pmass0=1.36;
LiInv=(Slipp-Slipp0)+250.*(1-(PmassN./Pmass0));
LiInv=[0,LiInv];
Slipp=[Slipp0,Slipp];
PmassN=[Pmass0,PmassN];

time1=[0,time];

if mode0==1
saveas(gcf,fullfile(fpath,sprintf('%d_matrix.jpg',CellID)));
else
end


if mode0==0 || mode0==2
figure(3)
   scatter(time1,Slipp, 'o', 'MarkerFaceColor', 'b')
    xlim([0 22000])
    ylim([0 60])
    xlabel('Time (h)')
    ylabel('Slippage (mAh)')
    ax = gca;
    ax.XRuler.Exponent = 0;
    box on
    saveas(gcf,fullfile(fpath,sprintf('%d_slippage.jpg',CellID)));
figure(4)
   scatter(time1,PmassN, 'o', 'MarkerFaceColor', 'b')
    xlim([0 22000])
    ylim([0.9 1.5])
    xlabel('Time (h)')
    ylabel('Positive mass (g)')
    ax = gca;
    ax.XRuler.Exponent = 0;
      box on
    saveas(gcf,fullfile(fpath,sprintf('%d_PosMas.jpg',CellID)));
figure(5)
   scatter(time,Capacity, 'o', 'MarkerFaceColor', 'b')
    xlim([0 22000])
    ylim([150 250])
    xlabel('Time (h)')
    ylabel('Charge capacity (mAh)')
    ax = gca;
    ax.XRuler.Exponent = 0;
     box on
    saveas(gcf,fullfile(fpath,sprintf('%d_Cap.jpg',CellID)));  
figure(6)
   scatter(time,-NslippS, 'o', 'MarkerFaceColor', 'b')
    xlim([0 22000])
    ylim([0 10])
    xlabel('Time (h)')
    ylabel('Negative slippage (mAh)')
    ax = gca;
    ax.XRuler.Exponent = 0;
     box on
    saveas(gcf,fullfile(fpath,sprintf('%d_Nslippage.jpg',CellID)));
figure(7)
   scatter(time,NmassO, 'o', 'MarkerFaceColor', 'b')
    xlim([0 22000])
    ylim([0.95 1.1])
    xlabel('Time (h)')
    ylabel('Negative mass (g)')
    ax = gca;
    ax.XRuler.Exponent = 0;
     box on
saveas(gcf,fullfile(fpath,sprintf('%d_NMass.jpg',CellID))); 
figure(8)
   scatter(time1,LiInv, 'o', 'MarkerFaceColor', 'b')
    xlim([0 22000])
    ylim([0 60])
    xlabel('Time (h)')
    ylabel('Lithium inventory loss (mAh)')
    ax = gca;
    ax.XRuler.Exponent = 0;
     box on
saveas(gcf,fullfile(fpath,sprintf('%d_LiInv.jpg',CellID)));  
output=[time1',PmassN',Slipp',LiInv'];
csvwrite(fullfile(fpath,sprintf('%d_output1.csv',CellID)),output)
output=[time',-PslippM',-NslippS',NmassO',Capacity'];
csvwrite(fullfile(fpath,sprintf('%d_output2.csv',CellID)),output)



else
end



    
    
  %{
PXX=PXX';
PYY=PYY';
GG=GG';

PX=PXX(:);
PY=PYY(:);
GG0=GG(:);

GGG=[PX,PY,GG0];
    

   %csvwrite('GG.csv',GG)
   csvwrite('GGG.csv',GGG)
%}


   
