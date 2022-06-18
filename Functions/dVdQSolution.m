function [X,dV1,Xe,dVne] = dVdQSolution(DD,A0,B0,s,Prec1,Prec2,MM,NN,OO,SS)

M=MM;
N=NN;
O=OO;
S=SS;
%-------------------------------------------  
%EXPERIMENTAL


[H3,X,L3]=size(DD);


C1=DD(:,4,s); %Experimental Voltage
C2=DD(:,3,s); %Experimental Capacity

%------------------------------------------- 
%THEORETICAL

Nmass=0.98+O*0.02; %Negative active mass
Pmass=0.9+0.05*(N-1)/Prec1; %Positive active mass
Nslipp=-S; %Negative slippage
Pslipp=-2.5*M/Prec2; %Positive slippage



A(:,1)=A0(:,1);
A(:,2)=A0(:,2).*Nmass+Nslipp; %Negative capacity

[H1,L1]=size(A0);


B(:,1)=B0(:,1);
B(:,2)=B0(:,2).*Pmass+Pslipp; %Positive capacity

[H2,L2]=size(B0);

for j=1:H1-1
    for i=1:H2-1
   
        X(i)=B(i,2); %Positive Capacity
        P(i)=B(i,1); %Positive Voltage
        
        if (X(i)<=A(j+1,2) && X(i)>=A(j,2))
           
            V(i)=A(j,1)+(A(j+1,1)-A(j,1)).*((X(i)-A(j,2))./(A(j+1,2)-A(j,2))); %Negative voltage
            
            D(i) = P(i)-V(i); %Predicted full cell voltage
            
        end    
        
    end
end

  for i=2:H2-1 %Finish with H2-1, because D finish with H2-1 in the last i loop.
      
      dV1(i)=(D(i)-D(i-1))./(X(i)-X(i-1));
  end
  


for k=2:H3
    
    dVe(k)=(C1(k)-C1(k-1))./(C2(k)-C2(k-1));
end

dVe=transpose(dVe);

 dVne=[];
 Xe=[];

for k=1:H3-1
    for i=1:H2-1
        
        if (X(i)<=C2(k+1) && X(i)>=C2(k))
                     
            dVne(i)=dVe(k)+(dVe(k+1)-dVe(k)).*((X(i)-C2(k))./(C2(k+1)-C2(k))); %Ajusted experimental dV/dQ
            Xe(i)=X(i);
        end    
    end
end

   
  
   
  

   
   %{ 
   [H5,L5]=size(dVne);
  xx = linspace(1,L5,L5);
   figure(1)
    
   plot(xx,dVne)
    xlim([8 250])
    ylim([0 0.1])
    xlabel('Capacity(mAh)')
    ylabel('dV/dQ')
   %}
   
 %{
    plot(X,dV1)
    xlim([8 250])
    hold on
   % plot(C2,dVe)
   plot(Xe,dVne)
    xlim([8 250])
    ylim([0 0.025])
    xlabel('Capacity(mAh)')
    ylabel('dV/dQ')
    hold off 
  %}
  
 
   %{
    figure(2)
    plot(Xe,G)
    xlim([8 250])
    xlabel('Capacity(mAh)')
    ylabel('Error')
    pause(0.1)
   %}
    
 %{
    [H5,L5]=size(B0(:,2));
    
    xx = linspace(1,H5,H5);
     figure(2)
    plot(xx,B0(:,2))
    xlim([8 250])
    ylim([-70 250])
    xlabel('Position')
    ylabel('X Capacity')
    pause(0.1)
  %}



end