%trackdata = matfile('trackdata.mat','Writable',false);
load('trackdata.mat');
load('Fahrzeugdaten.mat');
erg=matfile('GeschwindigkeitVorberechnung.mat','Writable',true);
Genauigkeit=0.1;
v_max=2*pi*n_max*r_dyn/Getriebeuebersetzung;
% r=abs(r);
% sigma = 50;
% sz = 300;    % length of gaussFilter vector
% x = linspace(-sz / 2, sz / 2, sz);
% gaussFilter = exp(-x .^ 2 / (2 * sigma ^ 2));
% gaussFilter = gaussFilter / sum (gaussFilter); % normalize
% 
% r = filter (gaussFilter,1, r);
%r= conv (r,gaussFilter,'same');
% v=0;
% a_y=Fzg_Gewicht*v^2/r(1,1);
v=zeros(Punkte,1);
%v_erg=zeros(Punkte,1);
a_y=zeros(Punkte,1);
Fy=zeros(Punkte,1);
Fz=zeros(Punkte,4); %hl,hr,vl,vr
F_radxy=zeros(Punkte,4);
Fzx=zeros(Punkte,4);
Fzy=zeros(Punkte,4);
Fzz=zeros(Punkte,4);
parfor i=1:Punkte
    diff=1;
    a_y_temp=0;
    v_zuvor=0;
    while diff>Genauigkeit
        %a_y=Fzg_Gewicht*v^2/r(i,1);
        Fzy(i,:)=[-a_y_temp*Fzg_Gewicht*COG_Hoehe/(Fzg_Spurweite*2),a_y_temp*Fzg_Gewicht*COG_Hoehe/(Fzg_Spurweite*2),-a_y_temp*Fzg_Gewicht*COG_Hoehe/(Fzg_Spurweite*2),a_y_temp*Fzg_Gewicht*COG_Hoehe/(Fzg_Spurweite*2)];
        Fzz(i,:)=[g*Fzg_Gewicht*(Fzg_Radstand-COG_Abstand_nach_hinten)*(Fzg_Spurweite-COG_Abstand_nach_links)/(Fzg_Radstand*Fzg_Spurweite),g*Fzg_Gewicht*(Fzg_Radstand-COG_Abstand_nach_hinten)*COG_Abstand_nach_links/(Fzg_Radstand*Fzg_Spurweite),g*Fzg_Gewicht*(COG_Abstand_nach_hinten)*(Fzg_Spurweite-COG_Abstand_nach_links)/(Fzg_Radstand*Fzg_Spurweite),g*Fzg_Gewicht*(COG_Abstand_nach_hinten)*COG_Abstand_nach_links/(Fzg_Radstand*Fzg_Spurweite)];
        Fz(i,:)=Fzx(i,:)+Fzy(i,:)+Fzz(i,:);
        F_radxy_temp=Reibwert_Rad*Fz(i,:);
        F_radxy(i,:)=F_radxy_temp;
        Fy(i,1)=F_radxy_temp(1,1)+F_radxy_temp(1,2)+F_radxy_temp(1,3)+F_radxy_temp(1,4);
        a_y_temp=Fy(i,1)/Fzg_Gewicht;
        a_y(i,1)=a_y_temp;
        v(i,1)=sqrt(abs(Fy(i,1)*r(i,1)/Fzg_Gewicht));
        diff=abs(v(i,1)-v_zuvor);
        v_zuvor=v(i,1);
        v(i,1)=min(v_zuvor,v_max);
    end
    a_y(i,1)=v(i,1)^2/r(i,1);
    Fzy(i,:)=[-a_y(i,1)*Fzg_Gewicht*COG_Hoehe/(Fzg_Spurweite*2),a_y(i,1)*Fzg_Gewicht*COG_Hoehe/(Fzg_Spurweite*2),-a_y(i,1)*Fzg_Gewicht*COG_Hoehe/(Fzg_Spurweite*2),a_y(i,1)*Fzg_Gewicht*COG_Hoehe/(Fzg_Spurweite*2)];
    Fz(i,:)=Fzx(i,:)+Fzy(i,:)+Fzz(i,:);
    F_radxy(i,:)=Reibwert_Rad*Fz(i,:);
    
end
plot (v,'y'), hold on, plot(r), hold off , axis([0 Punkte -0 40]);