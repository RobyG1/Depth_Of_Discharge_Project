clc;
clear;

Nmass=1.08; %Negative active mass
Pmass=1.23; % Positive active mass
Nslipp=-1; %Negative slippage
Pslipp=-32.5; %Positive slippage



%Input negative curve
A = readmatrix('NG.csv');
%A0 matrix is created from the Excel sheet

A(:,2)=A(:,2).*Nmass+Nslipp;

[H1,L1]=size(A);

%Input positive curve
B = readmatrix('622A.csv');

B(:,2)=B(:,2).*Pmass+Pslipp;

[H2,L2]=size(B);



for j=1:H1-1
    for i=1:H2-1
   
        X(i)=B(i,2);
        P(i)=B(i,1);
        
        if (X(i)<=A(j+1,2) && X(i)>=A(j,2))
           
            V(i)=A(j,1)+(A(j+1,1)-A(j,1)).*((X(i)-A(j,2))./(A(j+1,2)-A(j,2)));
            
            D(i) = P(i)-V(i);
            
        end
        
        
    end
end


 figure(1)
    plot(X,P)
    xlabel('Capacity (mAh)')
    ylabel('Voltage (V)')
    
    hold on
    
    plot(X,D)
    plot(X,V)
    
    hold off
    
    C=[X;D;V;P];
    C=transpose(C);
    
    csvwrite('FullCurve_66842.csv',C)
    
    Q1=0;
    Q2=0;
    
    %{
   for i=2:H2-2
       
       
      if abs(D(i)-3)<abs(D(i-1)-3)&&abs(D(i)-3)<abs(D(i+1)-3)
          Q1=X(i);
      else
      end
      
      if abs(D(i)-4.3)<abs(D(i-1)-4.3)&&abs(D(i)-4.3)<abs(D(i+1)-4.3)
          Q2=X(i);
      else
      end     
     
   end
   


%}