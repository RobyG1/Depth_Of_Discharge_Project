%{
clc;
clear;

%UL=1; %1 if upper cutoff cells, 0 if lower cutoff cells

%Load curves
A = readmatrix('66830_QC.txt');
B = readmatrix('66830_VQ.txt');
C = readmatrix('66830_IT.txt');
%}
function D = Findcheckup(UL,A,B,C)

Q=A(:,2); %Capacity
n=A(:,1); %Cycle number
V=B(:,2); %Voltage
I=[0;C(:,2)]; %Current
T=[0;C(:,1)]; %Time



%Find checkup cycles
if UL==1
    Cd=find(V<3.2 & abs(I)<25 & diff(abs([0;I]))==0);
    %
else
    Cd=find(V>3.95 & abs(I)<25 & diff(abs([0;I]))==0);
  
end

nc=n(Cd);

nck=unique(nc); %All the checkup cycle number
nck(diff(nck)==1) = []; %Get only one checkup cycle per checkup


D=[n,T,Q,V,I];

D=D(ismember(D(:,1),nck),:);
D = D(D(:,5)>0, :);
end


%{
figure(1)
scatter(D(:,3),D(:,4))
xlim([0.7e7 1e7])



%[L,H] = size(D);

%D = D(all(D,2),:); %Remove all rows with zeros


%}

