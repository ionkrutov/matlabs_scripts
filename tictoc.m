%Works; Elapsed time [between tic and toc] is 398.995829 seconds.

% addpath('/Users/redshift/Desktop/science/matlab codes/functions')

Waves.Begin=0;

[x,y]=meshgrid(0.0005:0.0005:0.75, 0.0005:0.0005:0.25);

interval=0.03/27.211;
Ee=interval:interval:8/27.211;
Ee=Ee';
capacity=length(Ee); 

l=0;
k=1/2;

lf=1;
kf=3/2;
J=2;
config='Jf12';
Jion=1/2;

lfin=0;
kfin=1/2;
Jfin=1;

Q=13;

znak=sign((-1)*(-1)^(0+1)*(-1)^(Jion+k)*w6j(l,1,1,1/2,k,Jion)); 

ImportEl='Initial_2p6_Final_2p5_J1_best.txt';

DiplEl=caseread(char(ImportEl));
info=DiplEl(1:1,: );
DiplEl=DiplEl(3:length(DiplEl),: );
DiplEl=str2num(DiplEl); 
DiplEl(:,all(DiplEl == 0))=[];
DiplElEn=DiplEl(:,2)/2; 

i=17;
fitReDiplEl = fit(DiplElEn,DiplEl(:,i),'nearestinterp');
fitImDiplEl = fit(DiplElEn,-DiplEl(:,i+1),'nearestinterp');
    
ImportDipEl="Intial_2s22p5_(" + config + 'l' + num2str(l) + ')_k' + num2str(2*k) + '2_Final_2s22p5_(' + config + 'l' + num2str(lf) + ')_k' + num2str(2*kf) + '2' + '_J' + num2str(J) + '.out'; %%there will be initial configuration in Jk coupling
        ContDiplEl=caseread(ImportDipEl);
        ContDiplEl=str2num(ContDiplEl);
        fitCont=fit( [x(:),y(:)], ContDiplEl(:), 'linearinterp');

ImportPhEl="phase" + config + 'l' + num2str(lf) + 'k' + num2str(2*kf) + '2.txt';        
        PhEl=caseread(ImportPhEl);
        PhEl=PhEl(2:length(PhEl),: );
        PhEl=str2num(PhEl);
        fitContPhase=fit(PhEl(:,1)/2,PhEl(:,2), 'linearinterp');

ImportDipEl2="Intial_2s22p5_(" + config + 'l' + num2str(lf) + ')_k' + num2str(2*kf) + '2_J' + num2str(J) + '_Final_2s22p5_(' + config + 'l' + num2str(lfin) + ')_k' + num2str(2*kfin) + '2' + '_J' + num2str(Jfin) + '.out';
        ContDiplEl2=caseread(ImportDipEl2);
        ContDiplEl2=str2num(ContDiplEl2);
        fitCont2=fit( [x(:),y(:)], ContDiplEl2(:), 'linearinterp');
        
 ImportPhEl2="phase" + config + 'l' + num2str(lfin) + 'k' + num2str(2*kfin) + '2.txt'; 
        PhEl2=caseread(ImportPhEl2);
        PhEl2=PhEl2(2:length(PhEl2),: );
        PhEl2=str2num(PhEl2);
        fitContPhase2=fit(PhEl2(:,1)/2,PhEl2(:,2), 'linearinterp');

lJk1=sqrt(2*l+1)*clebschgordan(l,0,1,0,lf,0)...
             *(-1)^(Jion-1/2+1+l)*sqrt((2*1+1)*(2*J+1)*(2*k+1)*(2*kf+1))*w6j(1/2,kf,J,1,1,k)*w6j(Jion,kf,lf,1,l,k);
lJk2=sqrt(2*lf+1)*clebschgordan(lf,0,1,0,lfin,0)...
     *(-1)^(Jion-1/2+J+lf)*sqrt((2*J+1)*(2*Jfin+1)*(2*kf+1)*(2*kfin+1))*w6j(1/2,kfin,Jfin,1,J,kf)*w6j(Jion,kfin,lfin,Jfin,lf,kf);

secondcontDownDown=zeros(capacity,1);
secondcont={0, 0, 0, secondcontDownDown}; 

Eout=interval:interval:8.5/27.211;
    clear polecos1
    clear polecos2
    polecos1(:,1)=cos(fitContPhase(Eout)+angle(igamma(lf+1-1j./sqrt(2*(Eout')),0))-angle(igamma(l+1-1j./sqrt(2*(Eout')),0))-angle(fitReDiplEl(Eout)+1j*fitImDiplEl(Eout)));
    polecos2(:,1)=cos(fitContPhase2(Eout)+angle(igamma(lfin+1-1j./sqrt(2*(Eout')),0))-angle(igamma(lf+1-1j./sqrt(2*(Eout')),0))-fitContPhase(Eout));

TimeCont2={0, 0, 0, TimeContQ4mwmw};
energymin=[0 0 0 40*interval];
energymax=[0 0 0 166*interval];
    
index=4;
    
E1=energymin(index):interval:energymax(index);
E2=energymin(index):interval:energymax(index);
[E1,E2]=meshgrid(E1,E2); % E1 rows, E2 cols
Cont1=fitCont(E1,E2);
Cont2=fitCont2(E1,E2);   

tic
    for index=4
    %parfor index=1:4
    Z=3;
    Zf=(energymin(index)+2*interval)/interval;
    for Ef=energymin(index)+interval:interval:energymax(index)-2*interval

        Y=2;
        
        Fzzz=TimeCont2{index}(Z,Z,Z)*abs(fitReDiplEl(Ef)+1j*fitImDiplEl(Ef));
        
        firstcontpv=0;
        X=1;
        for En=energymin(index):interval:energymax(index)
            if abs(Ef-En)>10^(-10)
                firstcontpv=firstcontpv+interval*(TimeCont2{index}(X,Z,Z)*abs(fitReDiplEl(En)+1j*fitImDiplEl(En))*Cont1(Z,X)-Cont1(Z,Z)*Fzzz)*(1/(En-Ef));
            else
                firstcontpv=firstcontpv+Fzzz*(Cont1(Z,Z)*log((energymax(index)-Ef)/(Ef-energymin(index)))+lJk1*sign(l-lf)*sqrt(2*Ef)*polecos1(Z,1));                   
            end
        X=X+1;
        end
                
        for Ek=energymin(index)+interval:interval:energymax(index)-interval
            
            X=1;
                       
            Fyyz=TimeCont2{index}(Y,Y,Z)*abs(fitReDiplEl(Ek)+1j*fitImDiplEl(Ek));
            
            firstcont=0; 
            for En=energymin(index):interval:energymax(index)
                                  
                if abs(Ek-En)>10^(-10)
                    firstcont=firstcont+interval*(TimeCont2{index}(X,Y,Z)*abs(fitReDiplEl(En)+1j*fitImDiplEl(En))*Cont1(Y,X)-Cont1(Y,Y)*Fyyz)*(1/(En-Ek));
                else
                    firstcont=firstcont+Fyyz*(Cont1(Y,Y)*log((energymax(index)-Ek)/(Ek-energymin(index)))+lJk1*sign(l-lf)*sqrt(2*Ek)*polecos1(Y,1));                   
                end
                
                X=X+1;
            end
             
                if abs(Ek-Ef)>10^(-10)
                    secondcont{index}(Zf,1)=secondcont{index}(Zf,1)+interval*(firstcont*Cont2(Z,Y)-firstcontpv*Cont2(Z,Z))*(1/(Ek-Ef));
                else
                    secondcont{index}(Zf,1)=secondcont{index}(Zf,1)+firstcontpv*(Cont2(Z,Z)*log((energymax(index)-interval-Ef)/(Ef-energymin(index)-interval))+lJk2*sign(lf-lfin)*sqrt(2*Ef)*polecos2(Z,1));
                end
                
            Y=Y+1;
        end
        Z=Z+1;
        Zf=Zf+1;
    end
    SecondContDipElTimeInt{index}(:,1)=znak*1/sqrt((2*Jfin+1)*(2*J+1)*3)*clebschgordan(1,0,1,0,J,0)*clebschgordan(1,0,J,0,Jfin,0)*secondcont{index}(:,1).*exp(1j*fitContPhase2(Ee))*(1j^(-lfin)).*exp(1j*angle(igamma(lfin+1-1j./sqrt(2*(Ee)),0)));
    end
    
WaveDownDown=['w' num2str(Q+4) '_' 'downdown' '_' config 'l' num2str(lfin) '_k' num2str(2*kfin) '2_J' num2str(Jfin)];

    if isfield(Waves,WaveDownDown)
    Waves.(WaveDownDown)=Waves.(WaveDownDown)+SecondContDipElTimeInt{4}(:,1);
    else
    Waves.(WaveDownDown)=SecondContDipElTimeInt{4}(:,1);  
    end    
toc

