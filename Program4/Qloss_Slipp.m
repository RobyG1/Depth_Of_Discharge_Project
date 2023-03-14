clc;
clear;

Nmass=0.974; %Negative active mass
Pmass=1.296; % Positive active mass

Ucut=4.1; %Upper cutoff
Lcut=3.77; % Lower cutoff

S0=12.4; %Initial slippage
Smax=100; % Maximum slippage
Prec=0.1; %Precision
K=Smax./Prec;

%Input negative curve
A0 = readmatrix('NG.csv');
%A0 matrix is created from the Excel sheet

A0(:,2)=A0(:,2).*Nmass;

[H1,L1]=size(A0);

%Input positive curve
B = readmatrix('622A.csv');

B(:,2)=B(:,2).*Pmass;

[H2,L2]=size(B);



for k=1:K
X=[];
V=[];
P=[];
A=A0;
S(k)=Prec*(k)+S0;

A(:,2)=A(:,2)+S(k);



for j=1:H1-1
    for i=1:H2-1
   
        X(i)=B(i,2);
        P(i)=B(i,1);
        
        if (X(i)<A(j+1,2) && X(i)>A(j,2))
           
            V(i)=A(j,1)+(A(j+1,1)-A(j,1)).*((X(i)-A(j,2))./(A(j+1,2)-A(j,2)));
            
            D(i) = P(i)-V(i);
            
        end    
    end
end

    
    Q1=0;
    Q2=0;
    
   for h=2:H2-2
       
       
      if Lcut>D(h)&&Lcut<D(h+1)
          Q1=X(h)+(X(h+1)-X(h)).*((Lcut-D(h))./(D(h+1)-D(h)));
      else
      end
      
      if Ucut>D(h)&&Ucut<D(h+1)
          Q2=X(h)+(X(h+1)-X(h)).*((Ucut-D(h))./(D(h+1)-D(h)));
      else
      end    
     
   end
   
   if k==1
       Y0=(Q2-Q1);
   else
   end
   
   Y(k)=Y0-(Q2-Q1);
   
end

figure(1)
    plot(S,Y)
    xlabel('Slippage (mAh)')
    ylabel('Capacity loss (mAh)')
    
    F=[S;Y];
    
    F=transpose(F);
    
   csvwrite('Qloss_slipp.csv',F)