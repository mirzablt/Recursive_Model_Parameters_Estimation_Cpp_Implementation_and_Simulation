%Specifikacija radne tacke
plant_mdl = 'mpc_cstr_plant';
op = operspec(plant_mdl);

%Pocetna vrijednost koncentracije reaktanta na ulazu u reaktor
op.Inputs(1).u = 10;
op.Inputs(1).Known = true;

%Pocetna vrijednost temperature reaktanta na ulazu u reaktor
op.Inputs(2).u = 298.15;
op.Inputs(2).Known = true;

%Pocetna vrijednost temperature rashladne tecnosti
op.Inputs(3).u = 298.15;
op.Inputs(3).Known = true;

%Racunanje pocetnih vrijednosti
[op_point, op_report] = findop(plant_mdl,op);

%Nominalne vrijednosti za stanja x, izlaze y i ulaze u
x0 = [op_report.States(1).x;op_report.States(2).x];
y0 = [op_report.Outputs(1).y;op_report.Outputs(2).y];		
u0 = [op_report.Inputs(1).u;op_report.Inputs(2).u;op_report.Inputs(3).u];

%Linearizacija modela objekta upravljanja u tacki pocetnih vrijednosti
sys = linearize(plant_mdl, op_point);

%CAi i CA se ne upotrebljavaju u MPC upravljanju pa ih treba iskljuciti
sys = sys(1,2:3);

% U imlementaciji MPC kontrolera koristi se  vremenski diskretni model procesa,
%pa kontinualni model procesa treba diskretizovati
Ts = 0.5;
plant = c2d(sys,Ts);

%Podesavanje kontolera
%MPC podesavamo u pocetnoj radnoj tacki. Kada je MPC pokrenut u adaptivnom modu,
%model procesa se azurira u toku rada. Specifikacija tipova signala upotrijebljenih u MPC-u
plant.InputGroup.MeasuredDisturbances = 1;
plant.InputGroup.ManipulatedVariables = 2;
plant.OutputGroup.Measured = 1;
plant.InputName = {'Ti','Tc'};
plant.OutputName = {'T'};

%MPC kontroler podesavamo sa difaltnim horizontima predikcije  i upravljanja
%(PredictionHorizon = 10, %ControlHorizon=2)
mpcobj = mpc(plant);

%Postaviti  nominalne vrijednosti kontrolera
mpcobj.Model.Nominal = struct('X',x0,'U',u0(2:3),'Y',y0(1),'DX',[0 0]);

%Ulaz i izlaz objekta imaju razlicite redove velicina, pa treba postaviti faktore skaliranja
Uscale = [30 50];
Yscale = 50;
mpcobj.DV.ScaleFactor = Uscale(1);
mpcobj.MV.ScaleFactor = Uscale(2);
mpcobj.OV.ScaleFactor = Yscale;

%Zbog fizikalnih ogranicenja sistema hlađenja ,prirastaj vremenske promjene
%rashladne tecnosti  TC  je ogranicena na 2 stepena u minuti
mpcobj.MV.RateMin = -2;
mpcobj.MV.RateMax = 2;

%Bus objekat
load('PSAU_BUS.mat')


