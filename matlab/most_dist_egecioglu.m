% This function implements the method proposed by Egecioglu (information
% Processing Letters, 1989), for finding the most distant points in a set
% (diameter of a set)

function [indQ,indQhat,r]=most_dist_egecioglu(Phi)

[n,p]=size(Phi);
flag=ones(1,p); %flag(i)=0 means that the i-th point has been excluded from the data set

P=Phi(:,1);
indP=1;
[rho,i2]=max(sqrt(sum((Phi-P*ones(1,p)).^2))); %Finding the most distant point of Phi from P
Q=Phi(:,i2);
indQ=i2;

[r,j2]=max(sqrt(sum((Phi-Q*ones(1,p)).^2))); %Finding the most distant point of Phi from Q
Qhat=Phi(:,j2);
indQhat=j2;

e=0;
j=1;

while(j<=p)||(e==0)
    Q_old=Q;
    indQ_old=indQ;
    Qhat_old=Qhat;
    indQhat_old=indQhat;
    r_old=r;
    flag(indP)=0; %Excluding P, Q from the data set
    flag(i2)=0;
    if(sum(flag)<=1)  %If almost all points have been considered
        break
    else
        Phat=Q+(r/rho)*(P-Q);
        P=0.5*(Phat+Q);  %Updating P
        
        [rho,i2]=max(sqrt(sum((Phi-P*ones(1,p)).^2)));
        Q=Phi(:,i2);
        
        [r,j2]=max(sqrt(sum((Phi-Q*ones(1,p)).^2)));
        Qhat=Phi(:,j2);
        
        if(r<=r_old)
            r=r_old;
            Q=Q_old;
            Qhat=Qhat_old;
            indQ=indQ_old;
            indQhat=indQhat_old;
            break
        else
            j=j+1;
        end
    end
end