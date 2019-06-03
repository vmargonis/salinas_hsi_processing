% This function computes the proposed MST statistic for a given data set.
% It takes as parameters the data set (VECTORS IN COLUMNS) and the number of sampled points,
% "k_max" and returns (a) the id numbers of the most distant points and (b)
% the most distant points themselves

function [most_dist,repre]=most_dist_repre(Phi,k_max)

[n,p]=size(Phi);

[indQ,indQhat,r]=most_dist_egecioglu(Phi); %Determining the diameter of the data set

most_dist=[indQ; indQhat];

t1=dist(Phi',Phi(:,most_dist(1)))';
t2=dist(Phi',Phi(:,most_dist(2)))';
% t1=sqrt(sum((Phi-Phi(:,most_dist(1))*ones(1,p)).^2));
% t2=sqrt(sum((Phi-Phi(:,most_dist(2))*ones(1,p)).^2));
% if n==1
%     t1=sqrt(sum((Phi-Phi(:,most_dist(1))).^2));
%     t2=sqrt(sum((Phi-Phi(:,most_dist(2))).^2));
% end   
temp=[t1; t2];

%%%%%%%%%%%%%%%%%%%%%%%%%

% tot_mean_dist_tot=sqrt((Phi(:,most_dist(1))-Phi(:,most_dist(2)))'*(Phi(:,most_dist(1))-Phi(:,most_dist(2))));
tot_mean_dist_tot=dist((Phi(:,most_dist(1)))',Phi(:,most_dist(2)))';

for k=3:k_max
    % For initialization
    min_temp=min(temp);
    [q,ind]=max(min_temp);    %This criterion is equivalent to the next one (in comments)
    
    temp_mean_dist=0;
    most_dist=[most_dist; ind];

    for q=1:length(most_dist)
%         temp_mean_dist=temp_mean_dist+ sqrt( sum((Phi(:,ind)-Phi(:,most_dist(q))).^2) );
        temp_mean_dist=temp_mean_dist+ dist(Phi(:,ind)',Phi(:,most_dist(q)))';
    end

    %Computation of the total mean distance
    temp=0;
    len=length(most_dist);
    for w1=1:len
        for w2=w1+1:len
%            temp=temp+sqrt((Phi(:,most_dist(w1))-Phi(:,most_dist(w2)))'* (Phi(:,most_dist(w1))-Phi(:,most_dist(w2))));
           temp=temp+dist((Phi(:,most_dist(w1)))',Phi(:,most_dist(w2)))';
        end
    end
    temp=temp/(len*(len-1)/2);
    
    tot_mean_dist_tot=[tot_mean_dist_tot temp];
    
    temp_mean_dist=temp_mean_dist/(length(most_dist)-1);
    
    temp=[];
    for q=1:length(most_dist)
%         temp=[temp; sqrt(sum((Phi-Phi(:,most_dist(q))).^2)) ];     %  dista(most_dist(q),:)];
        temp=[temp; dist(Phi',Phi(:,most_dist(q)))' ];     %  dista(most_dist(q),:)];
    end
    
end

most_dist;
repre=Phi(:,most_dist);
end

