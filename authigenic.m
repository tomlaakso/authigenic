
function metric = tcolumn(phase)

%Run parameters
deox = 1;
nx = 9;          %# dissolved species
ny = 5;          %# solid species
nr = 13;         %# reactions
nplt = 100;      %number of times to plot
zmax = 1000;          %tinker (from 300)
dz = 10.^[-1:.1:2];  %tinker (fom 10^1)
z = 0+cumsum(dz);
z = [z z(end)+max(dz):max(dz):zmax];
z = -z;
dz = [dz repmat(max(dz),1,length(z)-length(dz))]';
dzx = repmat(dz,1,nx);
dzy = repmat(dz,1,ny);    
rho = 2;         %solid density (g/cm^3)
sur = 10;        %surface area (m^2/g)
g = .15;          %methanogenically-available fraction
zeta = 2;        %scale for lability factor (0 = off)
phi0 = 0.5;      %porosity
w = .2/3.1/10^7; %sedimentation rate (cm/s)
D = 2.1*10^-5;   %diffusivity (cm2/s)

%Convenience terms
phi = linspace(phi0,phi0,length(z)); 
D = repmat(D,1,nx); 
D = repmat(D,length(z),1);
ps = phi./(1-phi);
phix = repmat(phi',1,nx);
phiy = repmat(phi',1,ny);
 
%index definitions: dissolved
ix = 1; %o2
id = 2; %DIC
is = 3; %SO4
it = 4; %H2S
iy = 5; %Fe2+
ih = 6; %H2
im = 7; %CH4
ic = 8; %Ca
ia = 9; %Alk

%index definitions: solid
io = 1; %Corg
iz = 2; %FeOOH
ip = 3; %FeS2
il = 4; %CaCO3
ig = 5; %methanogenically-available Corg

%Reaction stoichiometries
%reaction by species
Rx = zeros(nr,nx);
Ry = zeros(nr,ny);
%Respiration (normalized to O2)
Rx(1,ix) = -1;Rx(1,id) = 1;Ry(1,io) = -1;                  %oxic resp
Rx(2,iy) = 4;Rx(2,id) = 1;Ry(2,iz) = -4;Ry(2,io) = -1;     %iron reduc
Rx(3,is) = -.5;Rx(3,it) = .5;Rx(3,id) = 1;Ry(3,io) = -1;   %sulfate reduc
Rx(4,id) = .5;Rx(4,im) = .5;Ry(4,io) = -1;                 %methanogenesis
%Re-oxidation (normalized to O2)
Rx(5,ix) = -1;Rx(5,it) = -.5;Rx(5,is) = .5;                %sulfide
Rx(6,ix) = -1;Rx(6,iy) = -4;Ry(6,iz) = 4;                  %ferrous
Rx(7,ix) = -1;Rx(7,im) = -.5;Rx(7,id) = .5;                %methane
Rx(8,ix) = -1;Rx(8,ih) = -2;                               %hydrogen
Rx(9,is) = -.5;Rx(9,im) = -.5;Rx(9,id) = .5;Rx(9,it) = .5; %AOM
%Precipitation & dissolution (normalized to solid)
Rx(10,iy) = -1;Rx(10,it) = -2;Rx(10,ih) = 1;Ry(10,ip) = 1; %pyrite
Rx(11,ic) = -1;Rx(11,id) = -1;Ry(11,il) = 1;               %calcite precip
Rx(12,:) = -Rx(11,:);Ry(12,:) = -Ry(11,:);                 %calcite dissol
Rx(13,it) = -2;Rx(13,ih) = .5;Ry(13,iz) = -1;Ry(13,ip) = 1;%pyrite 2 
%Methane corrections
Ry(1:4,ig) = Ry(1:4,io);

%Rate constants
k1 = .0009/3.1/10^7;         %Katsev & Crowe 2015
k2 = k1;
k3 = k1;
k4 = k1;
k5 = 3.5;                    %Millero GCA1989
k6 = 2;                      %Luther et al. Front Microbio 2011
k7 = 2*10^5/3.1/10^7;        %Reeburg 2007 Chem Rev
k8 = 2.3/3.1/10^7;          
k9 = 2.7*10^5/3.1/10^7;      %Reebug 1991 Deep Sea Res,
k10 = .1;                    %Rickard & Luther 2007
k11 = 50*10^12/(3.6*10^19);  %global mass balance
k11 = k11/3.1/10^7;
k12 = k11;
k13 = 10^-4;                 %Raiswell & Canfield 1998


a1 = .0009/3.1/10^7;         %Katsev & Crowe 2015
a2 = a1;
a3 = a1;
a4 = a1;
a5 = 3.5;                    %Millero GCA1989
a6 = 2;                      %Luther et al. Front Microbio 2011
a7 = 2*10^5/3.1/10^7;        %Reeburg 2007 Chem Rev
a8 = 2.3/3.1/10^7;          
a9 = 2.7*10^6/3.1/10^7;      %Reebug 1991 Deep Sea Res,
a10 = .1;                    %Rickard & Luther 2007
a11 = 50*10^12/(3.6*10^19);  %global mass balance
a11 = a11/3.1/10^7;
a12 = a11;
a13 = 10^-4;                 %Raiswell & Canfield 1998

%M-M Half-maxima
km = ones(nx,1);
km(ix) = 1*10^-6;
km(is) = 200*10^-6; 
km = km/1000; %M to mol/cm3 
km(iy) = 10^-4;

%Carbonate system parameters, in M
%Zeebe & Wolf-Gladrow Appendix A
T = 5+273;
S = 35;
I = 19.924*S/(1000-1.005*S);
kc1 = exp( 2.83655-2307.1266/T-1.5529413*log(T)...
    -(0.207608410+ 4.0484/T)*sqrt(S)...
    +0.0846834*S-0.00654208*S^(3/2)+log(1-0.001005*S));
kc2 = exp(-9.226508-3351.6106/T-.2005743*log(T)...
    -(0.106901773+23.9722/T)*sqrt(S)...
    +0.1130822*S-0.00846934*S^(3/2)+log(1-0.001005*S));
ksp = 10^(-171.9065-.077993*T+2839.319/T+71.595*log10(T)...
    +(-.77712+.0028426*T+178.34/T)*sqrt(S)...
    -.07711*S+.0041249*S^1.5);
Rpdb = 1.1237*10^(-2);


%Isotopic parameters
ea = -.1190*(T-273)+12.09; %Equilibrium epsilon, bicarbonate-carbacid
eb = -.0569*(T-273)+8.53;  %Equilibrium epsilon, carbonate-carbacid
a = ea/1000*Rpdb;          %(convenience term)
b = eb/1000*Rpdb;          %(convenience term)
eo = -25;                  %epsilon, organic-seawater
em = -45;                  %epsilon, methane-organic

%Initial conditions
%Dissolved species (M)
x = [3 .95*2000 1*10^3 10^-6 100*10^-6 10^-6 10^-6 10.3*10^3 2000]*10^-6; %M
x(2) = 2100*10^-6;x(end) = 2190*10^-6; %Gives Om=1 for [Ca] = 10.3mM
x = repmat(x,length(z),1);
x = x*10^-3; %mole/cm3
%Solid species (moles / solid cm3)
y = [.10*rho/(12) .040*rho/(56) 10^-6 0];
y = [y y(1)*g]; %methanogen-available fraction
y = repmat(y,length(z),1);
%Initial conditions with depth; altered for faster spin-up
y(:,2) = y(:,2).*linspace(1,10^-10,size(y,1))'; %eliminates massive oxidant crash at depth
y(:,end) = y(:,1)*g;
y(2:end,4) = 10^-15; %eliminates catastrophic dissolution / alk shifts
%Upper boundary conditions
x0 = x(1,:);   %mole/cm3
y0 = y(1,:)*w; %moles/(cm^2)/s
x0base = x0;
y0base = y0;

%Isotopic boundaries
d0d = 0;
r0d = (d0d/1000+1)*Rpdb/(1+(d0d/1000+1)*Rpdb);
d0o = d0d+eo;   r0o = (d0o/1000+1)*Rpdb/(1+(d0o/1000+1)*Rpdb);
d0m = d0d+eo+em;r0m = (d0m/1000+1)*Rpdb/(1+(d0m/1000+1)*Rpdb);
d0l = d0d;      r0l = (d0l/1000+1)*Rpdb/(1+(d0l/1000+1)*Rpdb);
rid = repmat(r0d,length(z),1);
rio = repmat(r0o,length(z),1);
rim = repmat(r0m,length(z),1);
ril = repmat(r0l,length(z),1);

rebootfail = 1;
if phase>=20
    %For reboot
    phase = phase-20;
    phaseA = phase;
    try
     load(['boot' num2str(phase+1) '.mat']);
     phase = phaseA;
     rebootfail = 0;
    end
end

if rebootfail

if phase == 0
    %Cold start
    spinup = 1;
    tope = 0;stope = 0;
    tmax = zmax/w/2;
elseif phase == 1
    %Calculate initial alkalinity profile
    load(['boot' num2str(phase) '.mat']);
    phase = 1;
    spinup = 0;
    tope = 0;stope = 0;
    %x(:,ic) = (x(:,id)/.75-x0base(ia)-2*(x(:,iy)-x0base(iy))+2*(x(:,is)-x0base(is)))/2+x0base(ic);
    newalk = x0base(ia)...
          +2*(x(:,ic)-x0base(ic))...
          +2*(x(:,iy)-x0base(iy))...
          -2*(x(:,is)-x0base(is));
    x(:,ia) = newalk;
    if sum(x(:,ia)/2>=x(:,id))>0
     x(:,id) = 2*x(:,ia);
    end
    tmax = zmax/w/2;
elseif phase == 2
    %Lock fluxes and reservoirs for isotope calculations
    load(['boot' num2str(phase) '.mat']);
    olx = x-dx*dt;
    oly = y-dy*dt;
    phase = 2;
    spinup = 0;
    tope = 1;stope = 1;
    tmax = 3*zmax/w;
elseif phase == 3
    %Run full chemistry + isotopes
    load(['boot' num2str(phase-1) '.mat']);
    phase = 3;
    spinup = 0;
    tope = 1;stope = 0;
    tmax = 2*zmax/w;
elseif phase>=10
    %For continuation
    phase = phase-10;
    load(['boot' num2str(phase+1) '.mat']);
end    

tau = 0;               %time elapsed
mplt = 0;              %plots elapsed
chr = autumn(nplt+1);  %color
cnt = 0;               %counters
ecnt = 0;elast = 0;    %counters
muer = 0;              %counters

end

while tau<tmax

if deox
 x(:,ix) = 10^-100;
end

if ~stope
%%%Advection/diffusion
dxdz = (x(1:end-1,:)-x(2:end,:))./dzx(2:end,:);
dxdz = [(x0-x(1,:))./dzx(1,:) ; dxdz];
xx = D.*phix.^2.*dxdz; xx = [xx ; 0*xx(end,:)];
d2xdz2 =  (xx(1:end-1,:)-xx(2:end,:))./dzx;

dydz = (y(1:end-1,:)-y(2:end,:))./dzy(2:end,:);
dydz = [(y0/w-y(1,:))./dzy(1,:) ; dydz]*w;
end

%%%Isotope advection/diffusion
if tope & ~stope
gradd = (rid(1:end-1).*x(1:end-1,id)-rid(2:end).*x(2:end,id))./dz(2:end);
gradd = [(r0d*x0(id)-rid(1)*x(1,id))./dz(1) ; gradd];
xx = D(:,1).*phi'.^2.*gradd; xx = [xx ; 0*xx(end,:)];
grad2d =  (xx(1:end-1,:)-xx(2:end,:))./dz;
gradm = (rim(1:end-1).*x(1:end-1,im)-rim(2:end).*x(2:end,im))./dz(2:end);
gradm = [(r0m*x0(im)-rim(1)*x(1,im))./dz(1) ; gradm];
xx = D(:,1).*phi'.^2.*gradm; xx = [xx ; 0*xx(end,:)];
grad2m =  (xx(1:end-1,:)-xx(2:end,:))./dz;

grado = (rio(1:end-1).*y(1:end-1,io)-rio(2:end).*y(2:end,io))./dz(2:end);
grado = [(r0o*y0(io)/w-rio(1)*y(1,io))./dz(1) ; grado]*w;
gradl = (ril(1:end-1).*y(1:end-1,il)-ril(2:end).*y(2:end,il))./dz(2:end);
gradl = [(r0l*y0(il)/w-ril(1)*y(1,il))./dz(1) ; gradl]*w;
elseif stope
gradd = (rid(1:end-1).*olx(1:end-1,id)-rid(2:end).*olx(2:end,id))./dz(2:end);
gradd = [(r0d*x0(id)-rid(1)*olx(1,id))./dz(1) ; gradd];
xx = D(:,1).*phi'.^2.*gradd; xx = [xx ; 0*xx(end,:)];
grad2d =  (xx(1:end-1,:)-xx(2:end,:))./dz;
gradm = (rim(1:end-1).*olx(1:end-1,im)-rim(2:end).*olx(2:end,im))./dz(2:end);
gradm = [(r0m*x0(im)-rim(1)*olx(1,im))./dz(1) ; gradm];
xx = D(:,1).*phi'.^2.*gradm; xx = [xx ; 0*xx(end,:)];
grad2m =  (xx(1:end-1,:)-xx(2:end,:))./dz;

grado = (rio(1:end-1).*oly(1:end-1,io)-rio(2:end).*oly(2:end,io))./dz(2:end);
grado = [(r0o*y0(io)/w-rio(1)*oly(1,io))./dz(1) ; grado]*w;
gradl = (ril(1:end-1).*oly(1:end-1,il)-ril(2:end).*oly(2:end,il))./dz(2:end);
gradl = [(r0l*y0(il)/w-ril(1)*oly(1,il))./dz(1) ; gradl]*w;
end

%Carbonate system
rcip = repmat(0,length(z),1);
 
if ~stope
 sqr = x(:,ia);
 lin = (x(:,ia)-x(:,id))*kc1;
 con = (x(:,ia)-2*x(:,id))*kc1*kc2;
 quadarg = sqrt(lin.^2-4*sqr.*con);
 H = (-[lin lin]+[-quadarg quadarg])./2./[sqr sqr];
 H = H(:,1).*(real(H(:,1))>0 & imag(H(:,1))==0)...
    +H(:,2).*(real(H(:,2))>0 & imag(H(:,2))==0);
 H2CO3 = H.^2/kc1./(H+2*kc2).*x(:,ia)*1000;
 HCO3 = H./(H+2*kc2).*x(:,ia)*1000;
 CO3 = kc2./(H+2*kc2).*x(:,ia)*1000;
 Om = x(:,ic)*1000.*CO3/ksp; 
 pH = -log10(H);
 cs = [H2CO3 HCO3 CO3];
 end
 
 if tope
 cub = -H2CO3*a*b;
 sqr = rid.*x(:,id)*1000*a*b + H2CO3*(a+b+2*a*b)...
       + HCO3*b*(1-a) + CO3*a*(1-b);
 lin = -rid.*x(:,id)*1000*(a+b+2*a*b) - H2CO3*(1+a)*(1+b) ...
       - HCO3*((1-a)*(1+b)-a*b) - CO3*((1+a)*(1-b)-a*b);
 con = rid.*x(:,id)*1000*(1+a)*(1+b) - HCO3*a*(1+b) - CO3*b*(1+a);
 sqr = sqr./cub;
 lin = lin./cub;
 con = con./cub;
r1 = [];
for i = 1:length(cub)
    bit = roots([1 sqr(i) lin(i) con(i)]);
    bit = bit([real(bit)>0 & real(bit)<1 & imag(bit)==0]);
    r1 = [r1;bit];
end
 d1 = 1000*(r1./(1-r1)/Rpdb-1);
 d2 = d1+ea;
 dcip = d2+0.9;
 rcip = (dcip/1000+1)*Rpdb./(1+(dcip/1000+1)*Rpdb);
 end

%Redox ladder
if ~stope
ax1 = 0;ay1 = 1; ax2 = 1;ay2 = 0;  
am = (ay2-ay1)/(ax2-ax1);
 tink = 10;
mk1 = am.*(tink*x(:,ix)/km(ix)-ax1)+ay1;    mk1(mk1<0) = 0;
mk2 = am.*(tink*y(:,iz)/km(iy)-ax1)+ay1;    mk2(mk2<0) = 0;
mk3 = am.*(tink*x(:,is)/km(is)-ax1)+ay1;    mk3(mk3<0) = 0;
mk2 = min([mk1';mk2'])';
mk3 = min([mk1';mk2';mk3'])';
lab = 1-y(:,io)/(y0(io)/w);
lab = 10.^(-zeta*lab);

%Reaction rates
% moles / volume / time
% depth by reaction #
     Rt = zeros(length(z),nr);
     Rt(:,1) = k1*y(:,io).*(x(:,ix))./(km(ix)+x(:,ix)).*lab;
     Rt(:,2) = k2*y(:,io).*(y(:,iz))./(km(iy)+y(:,iz)).*lab;
     Rt(:,3) = k3*y(:,io).*(x(:,is))./(km(is)+x(:,is)).*lab;
     Rt(:,4) = k4*y(:,ig);
     Rt(:,5) = k5*x(:,ix).*x(:,it);
     Rt(:,6) = k6*x(:,ix).*x(:,iy);
     Rt(:,7) = k7*x(:,ix).*x(:,im);
     Rt(:,8) = k8*x(:,ix).*x(:,ih);
     Rt(:,9) = k9*x(:,is).*x(:,im);
     Rt(:,10) = k10*x(:,it).*x(:,iy);
     Rt(:,11) = k11*(Om-1).*(Om>1);
     Rt(:,12) = -k12*y(:,il)*(44+12+16*3).*(Om-1).*(Om<=1); 
     Rt(:,13) = k13*x(:,it).*y(:,iz);
       
%Ladder adjustments
Rt(:,2) = Rt(:,2).*mk1;
Rt(:,3) = Rt(:,3).*mk2;
Rt(:,4) = Rt(:,4).*mk3;

%Spin-up
if spinup
 Rt(:,11) = 0;
 Rt(:,12) = 0;
end

%Final reaction rates
phir = repmat(phi',1,4);
Rt(:,1:4) = Rt(:,1:4).*(1-phir)./phir;
Rxx = Rt*Rx;                 %depth by species
Ryy = Rt*Ry.*phiy./(1-phiy); %depth by species
Rm = Rt(:,4)*Ry(4,io).*phi'./(1-phi'); %NOTE SIGN!!!
gprop = y(:,ig)./y(:,io);
gprop = 1*gprop; %1 is default (i.e.: is it more labile?)
gprop = gprop-(gprop-1).*(gprop>1);
Ryy(:,ig) = (Ryy(:,io)-Rm).*(gprop)+Rm;

%Isotope reaction rates
if tope
ps = phi'./(1-phi');
dio = (rio./(1-rio)/Rpdb-1)*1000;
rgen = ((dio+em)/1000+1)*Rpdb./(1+((dio+em)/1000+1)*Rpdb);
rres = ((dio-em)/1000+1)*Rpdb./(1+((dio-em)/1000+1)*Rpdb);
Tid = rid.*x(:,id);
Til = ril.*y(:,il);
Tio = rio.*y(:,io); 
Tim = rim.*x(:,im);
delid = (-Rt(:,11).*rcip + Rt(:,12).*ril...
             + sum(Rt(:,1:3)')'.*rio...
             + .5*Rt(:,4).*rres...
             + .5*(Rt(:,7)+Rt(:,9)).*rim...
             + grad2d);
delil = (Rt(:,11).*rcip.*ps - Rt(:,12).*ril.*ps + gradl);
delio = (-(sum(Rt(:,1:4)')').*rio.*ps + grado);
delim = (.5*Rt(:,4).*rgen - .5*(Rt(:,7)+Rt(:,9)).*rim + grad2m);
end

%%% Eulerian update %%%
dx = (d2xdz2 + Rxx);     %depth by species
dy = (dydz + Ryy);       %depth by species

if deox
 dx(:,ix) = 0;
end

%Alkalinity uptdate
if ~spinup
 dx(:,ia) = (2*dx(:,ic)-2*dx(:,is)+2*dx(:,iy))*(1-spinup);
 q = x(:,id)-x(:,ia)/2;
 dq = (x(:,id)+dx(:,id)-(x(:,ia)+dx(:,ia))/2)-(x(:,id)-x(:,ia)/2);
end

%Time stepping
%Time step is dynamic. It is set to force the largest possible
%change in the prognostic species to be 1%. 
eul = .01^(1+ecnt/10);
tolz = 10^-50;
tolx = abs(x)>tolz;
toly = abs(y)>tolz;
dt1 = eul/max(max(abs(dx(tolx)./x(tolx))));
dt2 = eul/max(max(abs(dy(toly)./y(toly))));
dt3 = dt2;
if ~spinup
 dt3 = eul/max(max(abs(dq./q)));
end
dt4 = dt3;
if tope
 irats = [delid./Tid delil./Til delio./Tio delim./Tim];
 irats(isnan(irats) | irats==Inf) = 0;
 dt4 = eul^1/max(max(abs(irats)));
end
dtset = [dt1 dt2 dt3 dt4];
dt = min(dtset);

%Update
x = x+dx*dt;                %depth by species
y = y+dy*dt;                %depth by species
if tope
rid = (Tid+delid*dt)./x(:,id);%rid(isnan(rid)) = 0;
ril = (Til+delil*dt)./y(:,il);ril(y(:,il)==0) = 0;
rio = (Tio+delio*dt)./y(:,io);rio(y(:,io)==0) = 0;
rim = (Tim+delim*dt)./x(:,im);%rim(isnan(rim)) = 0;

did = (rid./(1-rid)/Rpdb-1)*1000;
dil = (ril./(1-ril)/Rpdb-1)*1000;
end

end
tau = tau+dt;

%%%Isotope reactions
if stope
ps = phi'./(1-phi');
dio = (rio./(1-rio)/Rpdb-1)*1000;
rgen = ((dio+em)/1000+1)*Rpdb./(1+((dio+em)/1000+1)*Rpdb);
rres = ((dio-em)/1000+1)*Rpdb./(1+((dio-em)/1000+1)*Rpdb);
Tid = rid.*olx(:,id);
Til = ril.*oly(:,il);
Tio = rio.*oly(:,io);
Tim = rim.*olx(:,im);
delid = (-Rt(:,11).*rcip + Rt(:,12).*ril...
             + sum(Rt(:,1:3)')'.*rio...
             + .5*Rt(:,4).*rres...
             + .5*(Rt(:,7)+Rt(:,9)).*rim...
             + grad2d);
delil = (Rt(:,11).*rcip.*ps - Rt(:,12).*ril.*ps + gradl);
delio = (-(sum(Rt(:,1:4)')').*rio.*ps + grado);
delim = (.5*Rt(:,4).*rgen - .5*(Rt(:,7)+Rt(:,9)).*rim + grad2m);
dt0 = eul^1/max(max(abs([delid./Tid delil./Til delio./Tio delim./Tim])));
rid = (Tid+delid*dt0)./(olx(:,id)+dx(:,id)*dt0);
ril = (Til+delil*dt0)./(oly(:,il)+dy(:,il)*dt0);
rio = (Tio+delio*dt0)./(oly(:,io)+dy(:,io)*dt0);
rim = (Tim+delim*dt0)./(olx(:,im)+dx(:,im)*dt0);
tau = tau-dt+dt0;
did = (rid./(1-rid)/Rpdb-1)*1000;
dil = (ril./(1-ril)/Rpdb-1)*1000;
end

try
if muer
    tmu = tmu+dt;
    Rmu = Rmu+Rt*dt/tmu;
end
end

%
if tau>mplt*tmax/nplt;
mplt = mplt+1;
disp(mplt-1);
clear phaseA;
if phase == 0
    save(['boot1.' num2str(zzz) '.mat']);
elseif phase == 1
    save(['boot2.' num2str(zzz) '.mat']);
elseif phase == 2
    save(['boot3.' num2str(zzz) '.mat']);
elseif phase == 3
    save(['boot4.' num2str(zzz) '.mat']);
end
end

if (mplt == nplt-1) & (muer == 0)
 muer = 1;tmu = 0;Rmu = Rt*0;   
else
 muer = 0;
end

cnt = cnt+1;
end

metric = .1/dt;

end
