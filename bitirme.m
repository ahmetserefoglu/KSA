close all;
clear;
%saha boyutlarý metre cinsinden x ve y koordinatlarý
xm=100;
ym=100;
%SÝnk'in x ve y koordinatlarý
sink.x=1.5*xm;
sink.y=0.5*ym;
%Alandaki düðümlerin sayisi 
n=200;
%Bir düðümde Optimal Seçim Olasýlðý ile Küme baþý olmasý
p=0.2;
intermediate=1;
% Enerji Modeli ( Jul tüm deðerler)
% Baþlangýç ??Enerji
Eo=0.1;
%Eelec=Etx=Erx
ETX=50*0.000000001;
ERX=50*0.000000001;
% Ýletim Amplifikatör tipleri
Efs=10*0.000000000001;
Emp=0.0013*0.000000000001;
% Veri Toplama Enerji
EDA=5*0.000000001;
%Hetereogeneity için  Geliþmiþ Deðerleri olandan yüzde düðüm
m=0.0;
%\alpha
a=1;
%maksimum döngü sayýsý
rmax=100;
%Computation of do
do=sqrt(Efs/Emp);
%Rastgele Sensör Aðý Oluþturma
figure(1);
hold off;
for i=1:1:n
    S(i).xd=rand(1,1)*xm;
    XR(i)=S(i).xd;
    S(i).yd=rand(1,1)*ym;
    YR(i)=S(i).yd;
    S(i).G=0;
% baþlangýçta hiçbir küme baþlarý orada sadece düðümleri vardýr
    S(i).type='N';
    temp_rnd0=i;
    %Normal Düðümlerin  Rastgele Seçimi
      if (temp_rnd0>=m*n+1) 
        S(i).E=Eo;
        S(i).ENERGY=0;
        plot(S(i).xd,S(i).yd,'o', 'MarkerSize', 5, 'MarkerFaceColor', 'g');
        hold on;
      end
    %Geliþmiþ Düðümler  Rastgele Seçim
      if (temp_rnd0<m*n+1)  
        S(i).E=Eo*(1+a);
        S(i).ENERGY=1;
        plot(S(i).xd,S(i).yd,'+', 'MarkerSize', 5, 'MarkerFaceColor', 'r');
        hold on;
      end
end
S(n+1).xd=sink.x;
S(n+1).yd=sink.y;
plot(S(n+1).xd,S(n+1).yd,'o', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
%ilk iterasyon
figure(1);
%cluster baþý sayacý
countCHs=0;
% Tur baþýna CHS için sayacý
rcountCHs=0;
cluster=1;
countCHs;
rcountCHs=rcountCHs+countCHs;
flag_first_dead=0;

for r=0:1:rmax
  if(mod(r, round(1/p) )==0)
    for i=1:1:n
        S(i).G=0;
        S(i).cl=0;
    end
  end

hold off;

%ölen düðümlerin sayýsý
dead=0;
%Ölü Geliþmiþ Düðümlerininsayýsý
dead_a=0;
%Ölü Normal Düðümlerininsayýsý
dead_n=0;
% Ýstasyonu Bazlar iletilen bit için sayaç ve Küme Heads
packets_TO_BS=0;
packets_TO_CH=0;
%Ýstasyonu Bazlar iletilen bit için sayaç ve Küme Heads
%Tur baþýna
PACKETS_TO_CH(r+1)=0;
PACKETS_TO_BS(r+1)=0;
figure(1);

for i=1:1:n
    %ölü düðüm varsa kontrol
    if (S(i).E<=0)
        plot(S(i).xd,S(i).yd,'red .', 'MarkerSize', 5, 'MarkerFaceColor', 'r');
        dead=dead+1;
        if(S(i).ENERGY==1)
            dead_a=dead_a+1;
        end
        if(S(i).ENERGY==0)
            dead_n=dead_n+1;
        end
        hold on;    
    end
    if S(i).E>0
        S(i).type='N';
        if (S(i).ENERGY==0)  
        plot(S(i).xd,S(i).yd,'o');
        end
        if (S(i).ENERGY==1)  
        plot(S(i).xd,S(i).yd,'+', 'MarkerSize', 5, 'MarkerFaceColor', 'g');
        end
        hold on;
    end
end
plot(S(n+1).xd,S(n+1).yd,'x');


STATISTICS(r+1).DEAD=dead;
DEAD(r+1)=dead;
DEAD_N(r+1)=dead_n;
DEAD_A(r+1)=dead_a;
plot(S(n+1).xd,S(n+1).yd,'o',  'MarkerSize', 5, 'MarkerFaceColor', 'r');
         plot(S(n+1).xd,S(n+1).yd,'x');  
%Ne zaman ilk düðümü ölür
if (dead==1)
    if(flag_first_dead==0)
        first_dead=r;
        flag_first_dead=1;
    end
end

countCHs=0;
cluster=1;
for i=1:1:n
   if(S(i).E>0)
   temp_rand=rand;     
   if ( (S(i).G)<=0)

 %Küme Baþkanlarý seçimi
 if(temp_rand<= (p/(1-p*mod(r,round(1/p)))))
            countCHs=countCHs+1;
            packets_TO_BS=packets_TO_BS+1;
            PACKETS_TO_BS(r+1)=packets_TO_BS;
            
            S(i).type='C';
            S(i).G=round(1/p)-1;
            C(cluster).xd=S(i).xd;
            C(cluster).yd=S(i).yd;
             plot(S(i).xd,S(i).yd,'k*',  'MarkerSize', 5, 'MarkerFaceColor', 'r');
            
            distance=sqrt( (S(i).xd-(S(n+1).xd) )^2 + (S(i).yd-(S(n+1).yd) )^2 );
            C(cluster).distance=distance;
            C(cluster).id=i;
            X(cluster)=S(i).xd;
            Y(cluster)=S(i).yd;
            cluster=cluster+1;
           %Enerji  hesaplanmasý daðýldý
           distance;
           if (distance>do)
               S(i).E=S(i).E- ( (ETX+EDA)*(4000) + Emp*4000*( distance*distance*distance*distance )); 
           end
            if (distance<=do)
               S(i).E=S(i).E- ( (ETX+EDA)*(4000)  + Efs*4000*( distance * distance )); 
            end
            Energy_disp(r+1) =  S(i).E;
        end     
    
    end
  end 
end

STATISTICS(r+1).CLUSTERHEADS=cluster-1;
CLUSTERHS(r+1)=cluster-1;
Ecalc=0;d1=0;d2=0;E1=0;E2=0;min=0;distance1=0;distance2=0;
%Normal Düðümler Ýliþkili Küme BaþkanýSeçim
for i=1:1:n
   if ( S(i).type=='N' && S(i).E>0 )
     if(cluster-1>=1)
       min_dis=sqrt( (S(i).xd-S(n+1).xd)^2 + (S(i).yd-S(n+1).yd)^2 );
       min_dis_cluster=1;
       for c=1:1:cluster-1
           temp=sqrt( (S(i).xd-C(c).xd)^2 + (S(i).yd-C(c).yd)^2 ) ;
           if ( temp<min_dis )
               min_dis=temp;
               min_dis_cluster=c;
           end
       end
       if(min_dis>do)
       Esmall=( ETX*(4000) + Emp*4000*( min_dis * min_dis * min_dis * min_dis)); 
       else
           Esmall=( ETX*(4000) + Efs*4000*( min_dis * min_dis)); 
       end
 for j=1:1:n
 if(i~=j && S(j).E>0)
 d1=sqrt( (S(i).xd-S(j).xd)^2 + (S(i).yd-S(j).yd)^2 ) ;
 d2=sqrt( (S(j).xd-C(min_dis_cluster).xd)^2 + (S(j).yd-C(min_dis_cluster).yd)^2 ) ;
 if(d1>do)
 E1=( ETX*(4000) + Emp*4000*( d1 * d1 * d1 * d1));
 else
 E1=( ETX*(4000) + Emp*4000*( d1 * d1));  
 end
 if(d2>do)
 E2=( ETX*(4000) + Emp*4000*( d2 * d2 * d2 * d2));
 else
 E2=( ETX*(4000) + Emp*4000*( d2 * d2));
 end
 Ecalc=E1+E2;
 if(Ecalc<min)
 min=Ecalc;
 intermediate=j;
 distance1=d1;
 distance2=d2;
 end
 end
 end
 if(Ecalc<Esmall)
 if(distance1>do)
        S(i).E=S(i).E- ( ETX*(4000) + Emp*4000*( distance1 * distance1 * distance1 * distance1));
        S(intermediate).E=S(intermediate).E-( (ERX + EDA)*4000 );
    else
        S(i).E=S(i).E- ( ETX*(4000) + Emp*4000*( distance1 * distance1));
        S(intermediate).E=S(intermediate).E-( (ERX + EDA)*4000 );
    end
    if(distance2>do)
        S(intermediate).E=S(intermediate).E-( ETX*(4000) + Emp*4000*( distance2 * distance2 * distance2 * distance2));
    else
        S(intermediate).E=S(intermediate).E-( ETX*(4000) + Emp*4000*( distance2 * distance2));
    end
       %Enerji iliþkili küme Baþkaný 
    %        min_dis;
     %       if (min_dis>do)
      %          S(i).E=S(i).E- ( ETX*(4000) + Emp*4000*( min_dis * min_dis * min_dis * min_dis)); 
       %     end
        %    if (min_dis<=do)
         %       S(i).E=S(i).E- ( ETX*(4000) + Efs*4000*( min_dis * min_dis)); 
          %  end
        %Energy dissipated
 else
if (min_dis>do)
                S(i).E=S(i).E- ( ETX*(4000) + Emp*4000*( min_dis * min_dis * min_dis * min_dis)); 
            end
            if (min_dis<=do)
                S(i).E=S(i).E- ( ETX*(4000) + Efs*4000*( min_dis * min_dis)); 
            end    
 end
        if(min_dis>0)
            distance=sqrt( (S(C(min_dis_cluster).id).xd-(S(n+1).xd) )^2 + (S(C(min_dis_cluster).id).yd-(S(n+1).yd) )^2 );
          S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E- ( (ERX + EDA)*4000 ); 
           if (distance>do)
                S(C(min_dis_cluster).id).E=S(C(min_dis_cluster).id).E- ( (ETX+EDA)*(4000) + Emp*4000*( distance*distance*distance*distance )); 
            end
            if (distance<=do)
                S(C(min_dis_cluster).id).E=S(C(min_dis_cluster).id).E- ( (ETX+EDA)*(4000)  + Efs*4000*( distance * distance )); 
            end
          PACKETS_TO_CH(r+1)=n-dead-cluster+1; 
        end

       S(i).min_dis=min_dis;
       S(i).min_dis_cluster=min_dis_cluster;
           
   end
 end
end
hold on;

countCHs;
rcountCHs=rcountCHs+countCHs;
sum=0;
for i=1:1:n
if(S(i).E>0)
    sum=sum+S(i).E;
end
end
avg=sum/n;
STATISTICS(r+1).AVG=avg;
sum;

  warning('OFF');
[vx,vy]=voronoi(X,Y);
plot(X,Y,'r*',vx,vy,'b-');
hold on;
voronoi(X,Y);
axis([0 xm 0 ym]);

end
figure(2);
for r=0:1:24
     ylabel('Her düðüm Ortalama Enerji');
    xlabel('Tur Sayýsý');
    plot([r r+1],[STATISTICS(r+1).AVG STATISTICS(r+2).AVG],'blue');
    hold on;
end
figure(3);
for r=0:1:49
    ylabel('Her düðüm Ortalama Enerji');
    xlabel('Tur Sayýsý');
    plot([r r+1],[STATISTICS(r+1).AVG STATISTICS(r+2).AVG],'blue');
    hold on;
end
figure(4);
for r=0:1:74
     ylabel('Her düðüm Ortalama Enerji');
    xlabel('Tur Sayýsý');
    plot([r r+1],[STATISTICS(r+1).AVG STATISTICS(r+2).AVG],'blue');
    hold on;
end
figure(5);
for r=0:1:99
    ylabel('Her düðüm Ortalama Enerji');
    xlabel('Tur Sayýsý');
    plot([r r+1],[STATISTICS(r+1).AVG STATISTICS(r+2).AVG],'blue');
    hold on;
end
figure(6);
for r=0:1:24
       ylabel('Ölen Düðümlerin Sayisi');
    xlabel('Tur Sayisi');
    plot([r r+1],[STATISTICS(r+1).DEAD STATISTICS(r+2).DEAD],'blue');
    hold on;
end
figure(7);
for r=0:1:49
         ylabel('Ölen Düðümlerin Sayisi');
    xlabel('Tur Sayisi');
    plot([r r+1],[STATISTICS(r+1).DEAD STATISTICS(r+2).DEAD],'blue');
    hold on;
end
figure(8);
for r=0:1:74
       ylabel('Ölen Düðümlerin Sayisi');
    xlabel('Tur Sayisi');
    plot([r r+1],[STATISTICS(r+1).DEAD STATISTICS(r+2).DEAD],'blue');
    hold on;
end
figure(9);
for r=0:1:99
        ylabel('Ölen Düðümlerin Sayisi');
    xlabel('Tur Sayisi');
    plot([r r+1],[STATISTICS(r+1).DEAD STATISTICS(r+2).DEAD],'blue');
    hold on;
end


%  DEAD : Ölü düðümlerin sayýsýnýn rmax x 1 dizi / yuvarlak							  %
%  DEAD_A : Ölü Geliþmiþ düðümler / turun sayýsýnýn rmax x 1 dizi				  %
%  DEAD_N : Ölü Normal düðümler / turun sayýsýnýn rmax x 1 dizi                    %
% CLUSTERHS : Küme Baþkanlarý sayýsýnýn rmax x 1 dizi / yuvarlak                     %
% PACKETS_TO_BS : sayý paketlerin bir rmax x 1 dizi Baz Ýstasyonu / tura göndermek    %
%  PACKETS_TO_CH : paketlerin sayýsýnýn rmax x 1 dizi ClusterHeads göndermek / yuvarlak   %
%  first_dead : Ýlk düðüm öldü yuvarlak
   