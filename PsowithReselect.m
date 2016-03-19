
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                      %
% This is the LEACH [1] code we have used.                             %
% The same code can be used for FAIR if m=1                            %
%                                                                      %
% [1] W.R.Heinzelman, A.P.Chandrakasan and H.Balakrishnan,             %
%     "An application-specific protocol architecture for wireless      % 
%      microsensor networks"                                           % 
%     IEEE Transactions on Wireless Communications, 1(4):660-670,2002  %
%                          www.forum.wsnlab.ir                         %                    
%     Homaei@wsnlab.ir &  Farhadi@wsnlab.ir & Ranjbaran@wsnlab.ir      %
%                                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Field Dimensions - x and y maximum (in meters)
xm=100;
ym=100;



%x and y Coordinates of the Sink
sink.x=0.5*xm;
sink.y=0.5*ym;

%Number of Nodes in the field
n=100

%Optimal Election Probability of a node
%to become cluster head
p=0.1;

%Energy Model (all values in Joules)
%Initial Energy 
Eo=0.5;
%Eelec=Etx=Erx
ETX=50*0.000000001;
ERX=50*0.000000001;
%Transmit Amplifier types
Efs=10*0.000000000001;
Emp=0.0013*0.000000000001;
%Data Aggregation Energy
EDA=5*0.000000001;
%my parameters
c=2;
winit=0.4;
tempo=double(5);
distsave=double(100);
positioner=size(100);
ad=0;
cluspos=size(5);
clusposprev=size(5);
%maximum number of rounds
rmax=2000;
Etot=Eo*100;
Total_Energy=Etot;
Eleft=Total_Energy;
 E_round_dissipate=0;
%%%%%%%%%%%%%%%%%%%%%%%%% END OF PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%

%Computation of do
do=sqrt(Efs/Emp);

%Creation of the random Sensor Network
figure(1);
for i=1:1:n
    S(i).xd=rand(1,1)*xm;
    XR(i)=S(i).xd;
    S(i).yd=rand(1,1)*ym;
    YR(i)=S(i).yd;
    S(i).G=0;
    %initially there are no cluster heads only nodes
    S(i).type='N';
    U(i).xd=S(i).xd;
    U(i).yd=S(i).yd;
    U(i).velocityx=0;
    U(i).velocityy=0;
    U(i).pbestx=U(i).xd;
    U(i).pbesty=U(i).yd;
    U(i).pbestcost=10000;
    U(i).gbestcost=U(i).pbestcost;
    U(i).cost=0;
    U(i).Sposition=i;
    U(i).gbestx=U(i).xd;
    U(i).gbesty=U(i).yd;
        S(i).E=Eo;
        U(i).E=S(i).E;
        S(i).ENERGY=0;
        plot(S(i).xd,S(i).yd,'o');
        hold on;
   
    
end
clusposprev=[0,0,0,0,0];
S(n+1).xd=sink.x;
S(n+1).yd=sink.y;
S(n+1).E=500;
plot(S(n+1).xd,S(n+1).yd,'x');
    
        
% Before First Iteration (Initialization phase)
figure(1);

%counter for CHs
countCHs=0;
%counter for CHs per round
rcountCHs=0;
cluster=1;
for i=1:1:5
    cpos=randi(100);
    countCHs=countCHs+1;
   % packets_TO_BS=packets_TO_BS+1;
    %PACKETS_TO_BS(r+1)=packets_TO_BS;
            
    S(cpos).type='C';
    C(cluster).xd=S(cpos).xd;
    C(cluster).yd=S(cpos).yd;
    plot(S(cpos).xd,S(cpos).yd,'k*');
    distance=sqrt( (S(cpos).xd-(S(n+1).xd) )^2 + (S(cpos).yd-(S(n+1).yd) )^2 );
    C(cluster).distance=distance;
    C(cluster).id=cpos;
    X(cluster)=S(cpos).xd;
    Y(cluster)=S(cpos).yd;
    cluster=cluster+1;
            
    %Calculation of Energy dissipated
    distance;
    if (distance>do)
        Edissipate=( (ETX+EDA)*(4000) + Emp*4000*( distance*distance*distance*distance ));
        S(cpos).E=S(cpos).E- ( (ETX+EDA)*(4000) + Emp*4000*( distance*distance*distance*distance )); 
    end
    if (distance<=do)
        Edissipate=( (ETX+EDA)*(4000) + Efs*4000*( distance*distance ));
        S(cpos).E=S(cpos).E- ( (ETX+EDA)*(4000)  + Efs*4000*( distance * distance )); 
    end
    E_round_dissipate=E_round_dissipate+Edissipate;
end
for i=1:1:n
   if ( S(i).type=='N' && S(i).E>0 )
       if(cluster-1>=0)
           min_dis=10000;
           for c=1:1:cluster-1
               temp=sqrt( (S(i).xd-C(c).xd)^2 + (S(i).yd-C(c).yd)^2 ) ;
               if ( temp<min_dis )
                   min_dis=temp;
                   min_dis_cluster=c;
               end
           end
           S(i).min_dist=min_dis;
           S(i).Cposition=min_dis_cluster;
           U(i).gbestx=C(min_dis_cluster).xd;
           U(i).gbesty=C(min_dis_cluster).yd;
       
           %    Energy dissipated by associated Cluster Head
           min_dis;
           if (min_dis>do)
               Edissipate=( (ETX+EDA)*(4000) + Emp*4000*( min_dis * min_dis * min_dis * min_dis ));
               S(i).E=S(i).E- ( (ETX+EDA)*(4000) + Emp*4000*( min_dis * min_dis * min_dis * min_dis)); 
                
           end
           if (min_dis<=do)
               Edissipate=((ETX+EDA)*(4000) + Efs*4000*( min_dis * min_dis ));
               S(i).E=S(i).E- ( (ETX+EDA)*(4000) + Efs*4000*( min_dis * min_dis)); 
               
           end
           E_round_dissipate=E_round_dissipate+Edissipate;
           %Energy dissipated
           if(min_dis>0)
               Edissipate  =( (ERX + EDA)*4000 );
               S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E- ( (ERX + EDA)*4000 );  
               E_round_dissipate=E_round_dissipate+Edissipate;
         
           end
       end
   end
end
countCHs;
rcountCHs=rcountCHs+countCHs;
flag_first_dead=0;
flag_last_dead=0;
flag_mid_dead=0;

for r=0:1:rmax
    r
E_round_dissipate=0;
  %Operation for epoch
 

hold off;

%Number of dead nodes
dead=0;
%Number of dead Advanced Nodes
dead_a=0;
%Number of dead Normal Nodes
dead_n=0;

%counter for bit transmitted to Bases Station and to Cluster Heads
packets_TO_BS=0;
packets_TO_CH=0;
%counter for bit transmitted to Bases Station and to Cluster Heads 
%per round
PACKETS_TO_CH(r+1)=0;
PACKETS_TO_BS(r+1)=0;

%Loop for velocity and position update
if(dead<=95)
for ids=1:1:n
    
        c1=c+rand(1,1);
        c2=c1;
        ra=rand(1,1);
        sdf=randi(50);
        w=winit+(sdf/100);
        U(ids).velocityx=(w*U(ids).velocityx)+(c1*ra*(U(ids).pbestx-U(ids).xd))+(c2*ra*(U(ids).gbestx-U(ids).xd));
        U(ids).velocityy=(w*U(ids).velocityy)+(c1*ra*(U(ids).pbesty-U(ids).yd))+(c2*ra*(U(ids).gbesty-U(ids).yd));
        U(ids).xd=U(ids).xd+U(ids).velocityx;
        U(ids).yd=U(ids).yd+U(ids).velocityy;
        
        if(U(ids).xd>100)
            U(ids).xd=100;
        end
        if(U(ids).yd>100)
            U(ids).yd=100;
        end
        if(U(ids).xd<0)
            U(ids).xd=0;
        end
        if(U(ids).yd<0)
            U(ids).yd=0;
        end
    
end
%cost calculation
    nsum=0;
    csum=0;
    for kl=1:1:n
        if(S(kl).E>0)
            if(S(kl).type=='N')
                nsum=nsum+S(kl).E;
            else
                csum=csum+S(kl).E;
            end
        end
    end
    f1=nsum/csum;
    
    %loop to find f2 for each node
    for c=1:1:5
        count=0;
        for l=1:1:n
            if ( U(l).E>0 )
                bn=U(l).Sposition;
                if(S(bn).type=='N')
                    if(S(bn).Cposition==c)
                    count=count+1;
                    end
                end
            end
        end
        C(c).counter=count;
    end
    for ml=1:1:n
        if(U(ml).E>0)
            mn=U(ml).Sposition;
            di=0;
            for c=1:1:5
                
                if ( S(mn).type=='N' && S(mn).E>0 )
                    if(S(mn).Cposition==c)
                    
                        di=di+S(mn).min_dist;
                    else
                        di=di+(sqrt( (S(mn).xd-C(c).xd)^2 + (S(mn).yd-C(c).yd)^2 ));
                    
                    end
                end
            end
            for cl=1:1:5
                tempo(cl)=di/C(cl).counter;
            end
            f2=max(tempo);
            %cost calculation for each node
            U(ml).cost=(0.6*f1)+(0.4*f2);
            if(U(ml).pbestcost>U(ml).cost)
                U(ml).pbestx=U(ml).xd;
                U(ml).pbesty=U(ml).yd;
                U(ml).pbestcost=U(ml).cost;
            end
            if(U(ml).gbestcost>U(ml).pbestcost)
                U(ml).gbestx=U(ml).pbestx;
                U(ml).gbesty=U(ml).pbesty;
                U(ml).gbestcost=U(ml).pbestcost;
            end
        end
    end
    for jl=1:1:n
        positioner(jl)=0;
    end
    for kk=1:1:n
        flag=0;
        for pk=1:1:n
            distsave(pk)=0;
        end
        for fg=1:1:n
            distsave(fg)=sqrt( (U(kk).xd-S(fg).xd)^2 + (U(kk).yd-S(fg).yd)^2 ) ;
        end
        while(flag==0)
            [ds,ad]=min(distsave);
            for alfa=1:1:n
                if(ad==positioner(alfa))
                    flag=0;
                    distsave(ad)=9999;
                    break;
                else
                    flag=1;
                end
            end
        end
        if(flag==1)
            positioner(kk)=ad;
            U(kk).Sposition=ad;
            U(kk).E=S(ad).E;
        end
        
    end
end
figure(1);

for i=1:1:n
    %checking if there is a dead node
    if (S(i).E<=0)
        plot(S(i).xd,S(i).yd,'red .');
        dead=dead+1;
       
        hold on;    
    end
    if (S(i).E>0)
        S(i).type='N';
        if (S(i).ENERGY==0)  
        plot(S(i).xd,S(i).yd,'o');
        end
        hold on;
    end
end
for er=1:1:5
    cluspos(er)=0;
end
plot(S(n+1).xd,S(n+1).yd,'x');
%continue here to find optimal clusterhead
    if(dead<=90)
        for qa=1:1:5
            temp=9999;
            clustercount=0;
            positer=0;
            for ws=1:1:n
                if(U(ws).E>0)
                    ups=U(ws).Sposition;
                    
                    if(S(ups).Cposition==qa)
                        clutercount=clustercount+1;
                        if(temp>U(ws).gbestcost)
                            temp=U(ws).gbestcost;
                            positer=ws;
                        end
                    end
                end
            end
            if(clustercount>0)
            cluspos(qa)=U(positer).Sposition;
            else
                gmd=0;
                while(gmd==0)
                    wer=randi(100);
                    if(S(wer).E>0)
                        cluspos(qa)=wer;
                        gmd=1;
                    else
                        gmd=0;
                    end
                end
            end
        end
        
    else
        if(dead<=95)
            coun=1;
            while(coun<=5)
                was=randi(100);
                if(S(was).E>0)
                    cluspos(coun)=was;
                    coun=coun+1;
                end
            end
        else
            cluspos(coun)=n+1;
            coun=coun+1;
        end
        
    end
    %reselect mechanism
    if(dead<90)
        
        if(clusposprev(1)==cluspos(1))
            firstcount=firstcount+1;
            if(firstcount==5)
                sad=cluspos(1);
                for jh=1:1:n
                    if(U(jh).Sposition==sad)
                        break;
                    end
                end
                U(jh).xd=S(sad).xd;
                U(jh).yd=S(sad).yd;
                U(jh).velocityx=0;
                U(jh).velocityy=0;
            end
        else
            firstcount=0;
        end
        
        if(clusposprev(2)==cluspos(2))
            secondcount=secondcount+1;
            if(secondcount==5)
                lad=cluspos(2);
                for lh=1:1:n
                    if(U(lh).Sposition==lad)
                        break;
                    end
                end
                U(lh).xd=S(lad).xd;
                U(lh).yd=S(lad).yd;
                U(lh).velocityx=0;
                U(lh).velocityy=0;
            end
        else
            secondcount=0;
        end
        
        if(clusposprev(3)==cluspos(3))
            thirdcount=thirdcount+1;
            if(thirdcount==5)
                mad=cluspos(3);
                for mh=1:1:n
                    if(U(mh).Sposition==mad)
                        break;
                    end
                end
                U(mh).xd=S(mad).xd;
                U(mh).yd=S(mad).yd;
                U(mh).velocityx=0;
                U(mh).velocityy=0;
            end
        else
            thirdcount=0;
        end
        
        if(clusposprev(4)==cluspos(4))
            fourthcount=fourthcount+1;
            if(fourthcount==5)
                bad=cluspos(4);
                for bh=1:1:n
                    if(U(bh).Sposition==bad)
                        break;
                    end
                end
                U(bh).xd=S(bad).xd;
                U(bh).yd=S(bad).yd;
                U(bh).velocityx=0;
                U(bh).velocityy=0;
            end
        else
            fourthcount=0;
        end
        
        if(clusposprev(5)==cluspos(5))
            fifthcount=fifthcount+1;
            if(fifthcount==5)
                rad=cluspos(5);
                for rh=1:1:n
                    if(U(rh).Sposition==rad)
                        break;
                    end
                end
                U(rh).xd=S(rad).xd;
                U(rh).yd=S(rad).yd;
                U(rh).velocityx=0;
                U(rh).velocityy=0;
            end
        else
            fifthcount=0;
        end
        clusposprev=cluspos;
    end
STATISTICS(r+1).DEAD=dead;
DEAD(r+1)=dead;
%DEAD_N(r+1)=dead_n;
%DEAD_A(r+1)=dead_a;

%When the first node dies
if (dead==1)
    if(flag_first_dead==0)
        first_dead=r
        flag_first_dead=1;
    end
end
%when half nodes are dead
if (dead==50)
    if(flag_mid_dead==0)
        half_dead=r
        flag_mid_dead=1;
    end
end
%when last node dies
if (dead==100)
    if(flag_last_dead==0)
        last_dead=r
        flag_last_dead=1;
    
    end
end

countCHs=0;
cluster=1;
countnode=0;
for i=1:1:5
    posb=cluspos(i);
    if(dead<=95)
   if(S(posb).E>0)   
   

 %Election of Cluster Heads
 
            countCHs=countCHs+1;
            packets_TO_BS=packets_TO_BS+1;
            PACKETS_TO_BS(r+1)=packets_TO_BS;
            
            S(posb).type='C';
            C(cluster).xd=S(posb).xd;
            C(cluster).yd=S(posb).yd;
            plot(S(posb).xd,S(posb).yd,'k*');
            
            distance=sqrt( (S(posb).xd-(S(n+1).xd) )^2 + (S(posb).yd-(S(n+1).yd) )^2 );
            C(cluster).distance=distance;
            C(cluster).id=posb;
            X(cluster)=S(posb).xd;
            Y(cluster)=S(posb).yd;
            cluster=cluster+1;
            
            %Calculation of Energy dissipated
            distance;
            if (distance>do)
                Edissipate=( (ETX+EDA)*(4000) + Emp*4000*( distance*distance*distance*distance ));
                S(posb).E=S(posb).E- ( (ETX+EDA)*(4000) + Emp*4000*( distance*distance*distance*distance )); 
            end
            if (distance<=do)
                Edissipate=( (ETX+EDA)*(4000) + Efs*4000*( distance*distance ));
                S(posb).E=S(posb).E- ( (ETX+EDA)*(4000)  + Efs*4000*( distance * distance )); 
            end
            E_round_dissipate=E_round_dissipate+Edissipate;
      
    
    
   end 
    end
end


STATISTICS(r+1).CLUSTERHEADS=cluster-1;
CLUSTERHS(r+1)=cluster-1;

%Election of Associated Cluster Head for Normal Nodes
for i=1:1:n
   if ( S(i).type=='N' && S(i).E>0 )
       if(dead<=95)
     if(cluster-1>=0)
       min_dis=9999;
       
       for c=1:1:cluster-1
           temp=sqrt( (S(i).xd-C(c).xd)^2 + (S(i).yd-C(c).yd)^2 ) ;
           if ( temp<min_dis )
               min_dis=temp;
               min_dis_cluster=c;
           end
       end
       
       %Energy dissipated by associated Cluster Head
            min_dis;
            if (min_dis>do)
                Edissipate=( (ETX+EDA)*(4000) + Emp*4000*( min_dis * min_dis * min_dis * min_dis ));
                S(i).E=S(i).E- ( ETX*(4000) + Emp*4000*( min_dis * min_dis * min_dis * min_dis)); 
                
            end
            if (min_dis<=do)
                Edissipate=((ETX+EDA)*(4000) + Efs*4000*( min_dis * min_dis ));
                S(i).E=S(i).E- ( ETX*(4000) + Efs*4000*( min_dis * min_dis)); 
               
            end
            E_round_dissipate=E_round_dissipate+Edissipate;
        %Energy dissipated
        if(min_dis>0)
            Edissipate  =( (ERX + EDA)*4000 );
          S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E- ( (ERX + EDA)*4000 ); 
         PACKETS_TO_CH(r+1)=n-dead-cluster+1; 
         E_round_dissipate=E_round_dissipate+Edissipate;
         
        end
S(i).E;
       S(i).min_dis=min_dis;
       S(i).min_dis_cluster=min_dis_cluster;
           
     end
       else
         min_dis=sqrt( (S(i).xd-S(n+1).xd)^2 + (S(i).yd-S(n+1).yd)^2 );
         
       %Energy dissipated by associated Cluster Head
            min_dis;
            if (min_dis>do)
                Edissipate=( (ETX+EDA)*(4000) + Emp*4000*( min_dis * min_dis * min_dis * min_dis ));
                S(i).E=S(i).E- ( ETX*(4000) + Emp*4000*( min_dis * min_dis * min_dis * min_dis)); 
                
            end
            if (min_dis<=do)
                Edissipate=((ETX+EDA)*(4000) + Efs*4000*( min_dis * min_dis ));
                S(i).E=S(i).E- ( ETX*(4000) + Efs*4000*( min_dis * min_dis)); 
               
            end
            E_round_dissipate=E_round_dissipate+Edissipate;
        %Energy dissipated
        if(min_dis>0)
            Edissipate  =( (ERX + EDA)*4000 );
          S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E- ( (ERX + EDA)*4000 ); 
         PACKETS_TO_CH(r+1)=n-dead-cluster+1; 
         E_round_dissipate=E_round_dissipate+Edissipate;
         
        end
         
       end
   end
end
hold on;

countCHs;
rcountCHs=rcountCHs+countCHs;
if(r==0)
        Total_Energy(r+1)=Total_Energy-E_round_dissipate;
else
        Total_Energy(r+1)=Total_Energy(r)-E_round_dissipate;
end
 if(r==0)
           PACKETS_TO_BS_TILL_ROUND(r+1)=PACKETS_TO_BS(r+1);
 else
           PACKETS_TO_BS_TILL_ROUND(r+1)=PACKETS_TO_BS(r+1)+PACKETS_TO_BS_TILL_ROUND(r);
 end
Eleft=Eleft-E_round_dissipate;


%Code for Voronoi Cells
%Unfortynately if there is a small
%number of cells, Matlab's voronoi
%procedure has some problems

%[vx,vy]=voronoi(X,Y);
%plot(X,Y,'r*',vx,vy,'b-');
% hold on;
% voronoi(X,Y);
% axis([0 xm 0 ym]);

end
first_dead

last_dead
for i=1:r+1
    round(i)=i;
    aliveleach(i)=n-STATISTICS(i).DEAD;
end
figure(5);
hold on;
plot(round,aliveleach,'g');
xlabel('Number of rounds --->');
ylabel('Number of nodes alive --->');
h=legend('leach','pso',2);
set(h,'Location','NorthEast');


figure(6);
hold on;
y=1:r+1;
plot(y,Total_Energy,'g');
xlabel('Number of Rounds --->');
ylabel('Total Energy of all the Clusters --->');
h=legend('leach','pso',2);
set(h,'Location','NorthEast');

hold on;
%plot(first_dead,Total_Energy(first_dead+1),'bo');

figure(7);
hold on;
plot(round,PACKETS_TO_BS_TILL_ROUND,'g');
hold on;
%plot(first_dead,PACKETS_TO_BS_TILL_ROUND(first_dead+1),'bo');
xlabel('Number of Rounds --->');
ylabel('Number of Packets Transmitted to Base Station --->');
h=legend('leach','pso',2);
set(h,'Location','NorthEast');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   STATISTICS    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                     %
%  DEAD  : a rmax x 1 array of number of dead nodes/round 							  %
%  DEAD_A : a rmax x 1 array of number of dead Advanced nodes/round					  %
%  DEAD_N : a rmax x 1 array of number of dead Normal nodes/round                     %
%  CLUSTERHS : a rmax x 1 array of number of Cluster Heads/round                      %
%  PACKETS_TO_BS : a rmax x 1 array of number packets send to Base Station/round      %
%  PACKETS_TO_CH : a rmax x 1 array of number of packets send to ClusterHeads/round   %
%  first_dead: the round where the first node died                                    %
%                                                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






