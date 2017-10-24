%Einlesen der x und y Koordinaten
trackdata = matfile('trackdata.mat','Writable',true);
xy=importdata('xy.txt');
x=xy(1,:);
y=xy(2,:);

%Erstellung des Splines
df = diff(xy,1,2);
t = cumsum([0, sqrt([1 1]*(df.*df))]);
track = csapi(t,xy);

%Erstellung der Spline-Ableitungen
dtrack=fnder(track); %erste Ableitung
ddtrack=fnder(dtrack); %zweite Ableitung

%Berechnung von Kappa und Radius an Abschnitten s
tracklength=round(t(length(t)),2); %hier 129.5021m Streckenlänge (letzte Spalte in Matrix t)
unterteilung=0.01; %in m
Punkte=round(tracklength/unterteilung,0);
s = [0:unterteilung:tracklength];
K=zeros(Punkte,1);
r=zeros(Punkte,1);
parfor i = 1:(tracklength/unterteilung)
    X=fnval(track,s(i));
    Y=fnval(dtrack,s(i));
    Z=fnval(ddtrack,s(i));
    K(i)=(Y(1,1)*Z(2,1)-Z(1,1)*Y(2,1))/(Y(1,1)^2+Y(2,1)^2)^(3/2);
    r(i)=(1/K(i));
end

%Plottausgabe vom Track/Kappa&Radius
%fnplt(dtrack), hold on, fnplt(ddtrack), fnplt(track), plot(x,y,'o'), hold off
%plot (K,'black'), hold on, plot (r,'.'), hold off, axis([0 Punkte -100 100]);

%Ausgabe in Datei trackdata.mat
trackdata.r=r; trackdata.k=K; trackdata.xy=xy; trackdata.track=track; trackdata.unterteilung=unterteilung; trackdata.tracklength=tracklength; trackdata.Punkte=Punkte;