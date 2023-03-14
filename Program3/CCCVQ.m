clc;
clear;
addpath('Input')
addpath('Functions')

CellID=81336;
HasCheckups=1; %Is 1 if the cell has checkup cycles and 0 if the cell has no checkup cycles
UpperorLower=0; %1 if fixed upper cutoff cells, 0 if fixed lower cutoff cells. 1 or 0 is fine for fixed middle cutoff cells. Keep the value at 1 if your cells are 100% DOD.
ScanNumber=9;
%maxV=4.2; %Maximum voltage in the data. Delete all possible noise above maxV.
%minV=2.9; %Minimum voltage in the data. Delete all possible noise below minV.



%Load curves
A = readmatrix(sprintf('%d_QC.txt',CellID)); %Input capacity vs cycle curves 
B = readmatrix(sprintf('%d_VQ.txt',CellID)); %Input capacity vs voltage curves
C = readmatrix(sprintf('%d_IT.txt',CellID)); %Input current vs time curves


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
%D=D(D(:,5)>0,:); %Only remove the discharge curves
end

%Clean data
%D = D(all(D,2),:); %Remove rows that have zeros
%D = D(D(:,4) <= maxV, :); %Remove all rows with V > maxV Volts
%D = D(D(:,4) >= minV, :); %Remove all rows with V < minV Volts

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

    
    
for mm=1:ScanNumber   
hold on

V=DD(:,4,mm);
Q=DD(:,3,mm);
n=DD(10,1,mm)

plot(Q,V)
ylim([3.9 4.2])

hold off
end
