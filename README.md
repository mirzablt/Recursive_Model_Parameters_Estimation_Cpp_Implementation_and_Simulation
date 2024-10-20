# Rekurzivne metode identifikacije sistema i njihova primjena u adaptivnim sistemima upravljanja

Ovaj repozitorij sadrži implementacije nekih od rekurzivnih algoritama parametarske identifikacije sistema koji su primjenjivi i u sistemima upravljanja koji u sebe uključuju identifikaciju objekta upravljanja. Algoritmi su implementirani u C++ programskom jeziku kao blok S-funkcije u Matlab Simulink-u i ovaj blok je korišten u simulacijama adaptivnog upravljanja.

## Uvod

Rekurzivne metode za estimaciju parametara sistema koriste se za identifikaciju matematičkih modela dinamičkih sistema na osnovu mjerenih ulazno-izlaznih podataka. Za razliku od klasičnih metoda gdje se svi podaci obrađuju istovremeno, rekurzivne metode ažuriraju parametre sistema u realnom vremenu dok pristižu novi podaci. Kod rekurzivnih metoda identifikacije, koje se još zovu i on-line metode identifikacije, estimati parametara se računaju rekurzivno u vremenu. Ovo znači da, ako postoji estimat koji je izračunat na osnovu izmjerenih podataka do prethodnog trenutka uzorkovanja, tada je estimat za tekući trenutak izračunat modifikacijom prethodnog estimata. Postoje različite varijante rekurzivnih algoritama identifikacije pogodnih za sisteme sa vremenski promjenjivim parametrima.

---

## Algoritmi estimacije parametara modela sistema

### Rekurzivni LS metod sa eksponencijalnim faktorom iščezavanja

Ovaj algoritam se koristi za online identifikaciju sistema i estimaciju parametara dinamičkih sistema koji se mijenjaju tokom vremena. Algoritam ažurira estimat parametara sistema na osnovu novopristiglih podataka i koristi faktor iščezavanja \( \lambda \) kako bi dao veću težinu novijim podacima, a smanjio utjecaj starijih.

### Rekurzivni LS metod sa konačnim vremenskim okvirom (Windowed RLS)

Ovaj algoritam koristi klizni prozor fiksne veličine i samo ograničen broj recentnih podataka za procjenu trenutnih parametara. Ažuriranje se vrši u dva koraka: "updating" i "downdating".

### Kalmanov filter

Kalmanov filter procjenjuje parametre modela kroz dvije faze: **predikciju** i **ažuriranje**. Kombinuje prethodne procjene sa novim mjerenjima, koristeći Kalmanovo pojačanje za optimizaciju procjena u prisustvu šuma.

### BFGS algoritam

BFGS metod je iterativna optimizacijska metoda koja se koristi za nelinearne sisteme kada se procjena parametara ne može riješiti u zatvorenoj formi. Funkcija cilja može biti suma kvadratnih grešaka između predviđenih i stvarnih vrijednosti izlaza sistema.

---

## Ulazi i izlazi algoritama

- **Ulazi**: Trenutne i prethodne vrijednosti ulaza i izlaza procesa, struktuirane kao regresori.
- **Izlazi**: Vektor parametara modela i greška predikcije. Ulazi i izlazi algoritma predstavljaju i ulaze i izlaze s-funkcije implementirane u Matlab Simulink-u.

---

## Linearni modeli

Ovi modeli opisuju odnose između ulaza i izlaza sistema pomoću diferencijalnih jednadžbi:
- **ARX**: Najjednostavniji model, koristi aditivni šum.
- **ARMAX**: Proširuje ARX dodavanjem pokretnog prosjeka za šum.
- **BJ**: Razdvaja dinamiku sistema i šuma.
- **OE**: Minimizira grešku bez eksplicitnog modeliranja šuma.

---

## Datoteke

- `ERLSLibrary.cpp`: U ovom fajlu su implementirane funkcije za struktuiranje izmjerenih podataka o sistemu u vektore (regresore) i matrice, te funkcije za dinamičku alokaciju memorije potrebne za smještanje izmjerenih podataka. Funkcije su detaljnije opisane u komentaru koda koji se nalazi u ovom fajlu;
- `BFGSAlgorithm.cpp`:U ovom fajlu je implementiran BFGS optimizacijski algoritam, funkcija cilja, algoritam linijskog pretraživanja i drugih funkcija koje su pozivane u navedenim implementacijama i koje su detaljno opisane u komentaru koda koji se nalazi u navedenom fajlu;
- `MatrixLibrary.cpp`: Operacije sa matricama i vektorima.
- `Recursive_Estimation_wrapper.cpp`: U ovom fajlu su implementirani gore opisani rekurzivni algoritmi, implementirna je logika koja omogučava izbor algoritma estimacije preko parametara bloka S-funkcije.
- `MPC_controller_parameters.m`: Pored učitavanja parametara MPC kontrolera korištenog u simulaciji, izvršavanjem ovog fajla učitavaju se inicijalne vrijednosti nekih parametara simulcije, što je detaljnije opisano u komentaru koda u ovom fajlu.
- `PSAU_BUS.mat`: Za proces sa vise ulaza, estimirani parametri sistema su struktuirani u matricu. Izlaz bloka S-funkcije u kojoj su implementirani algoritmi estimacije parametara sistema je sabirnica `PSAU_BUS`, koju čine parametri. U `PSAU_BUS.mat` su definisane dimenzije `A`, `B`, `F`, `C`, `D` parametara.
- `Documentation.pdf`: Ovaj projekat je rađen kao seminarski rad na Elektrotehničkom fakultetu u Sarajevu. U ovom fajlu je dokumentovan ovaj rad. U njemu su opisani detalji koji se tiču linearnih modela sistema, rekurzivne estimacije parametara modela i  primjene BFGS algoritma optimizacije u estimaciji parametara modela. Opisane su simulacije i rezultati simulacija. U komentaru koda se poziva na numerisane dijelove ovog dokumenta. Kod koji je dokumentovan u ovom fajlu je dodatno dorađen u smislu njegove preglednosti a komentar koda je opširniji i detaljniji.

---

## Simulacije

Da bi koristili ovaj projekt, potrebno je na računaru imati instaliran Matlab R2020a ili svježiju verziju. U simulacijama se koristi rekurzivna estimacija parametara modela objekta upravljanja (u zatvorenoj upravljačkoj konturi). Objekat upravljanja u primjerima je hemijski reaktor koji je izrazito nelinearan sistem. 
### Adaptivno MPC upravljanje sa real-time identifikacijom procesa
Simulacija Adaptivnog MPC (*Model Predictive Control*) upravljanja koristi ARX linearni model procesa. Kako je objekat nelinearan, to je linearni model objekta vremenski promjenjiv. Prema tome ovaj model  se mora  ažurirati u svakom trenutku odabiranja, za šta nam je služi  implemetirani blok rekurzivnog estimatora parametara, u kojem je (preko parametara bloka) odabran Kalmanov filter za metod estimacije parametara.
Da bi simulirali nelinearno MPC upravljanje i demonstrirali rad navedenog algoritama estimacije parametara, potrebno je:
1.	Instalirati C++ kompiler u Matlab-u; 
2.	Klonirati repozitorij;
3.	Otvoriti m-fajl: `MPC_controller_parameters.m` koji se nalazi u ovom repozitoriju. Pokrenuti izvršenje ovog fajla klikom na dugme “Run”  koje se nalazi u alatnoj traci Matlab-a. 
4.	Dvoklikom na `AdaptiveMPC_Simulation.slx` fajl otvara se Simulink blok dijagram koji predstavlja simulaciju sistema upravljanja sa on-line estimcijom parametara implementiranog u S-funkciji;
5.	 Dvoklikom na blok S-funkcije (nazvane: *Recursive Parameters Estimator*) i klikom na dugme *Build*, kompajlira se C++ kod implementiran u S-funkciji;
6.	Klikom na dugme *Run* koje se nalazi u alatnoj traci Simulink-a pokreće se simulacija sistema upravljanja.

### PI adaptivno upravljanje sa real-time identifikacijom procesa
PI adaptivno upravljanje u kojem  se estimacija parametara linearnog  modela procesa (objekta upravljanja) vrši samo pri promjeni radne tačke procesa. Zbog nelinearnosti objekta, za različite radne tačke imamo različite linearne modele, odnosno prenosne funkcije koje karakterišu estimirani parametri. Ovi parametri su upotrijebljeni za podešavanje pojačanja PI regulatora. U ovom primjeru je izabran ARMAX model, a metoda estimacije je ERLS sa faktorom iščezavanja. 
Da bi simulirali navedeno upravljanje i demonstrirali rad navedenog algoritama estimacije parametara, pored navedena prva dva koraka prethodnog primjera potrebno je:
1.	Dvoklikom na `PSAU_BUS.mat` fajl iz navedenog repozitorija, ovaj fajl se pojavljuje u radnom prostoru Matlab-a;
2.	Dvoklikom na `PSAU_BUS.mat` u radnom prostoru otvara se *Bus Editor* u kojem je vrijednost `B` sa `[2 3]` potrebno promijeniti na vrijednost `[1 3]`;
3.	Dvoklikom na `AdaptivePI_Simulation.slx` fajl otvara se Simulink blok dijagram koji predstavlja  simulaciju sistema upravljanja; 
4.	Dvoklikom implementirani blok `Recursive Online Estimator` (S-funkcija) otvara se prozor u kome je u mapi *S-function Parameters*, parametar bloka `nb`, sa `[2 2]` potrebno promijeniti na vrijednost `[2]`;
5.	Klikom na dugme *Build* u prozoru s-funkcije, kompajlira se C++ kod implementiran u s-funkciji.
6.	Klikom na dugme *Run* koje se nalazi u alatnoj traci Simulink-a u simulaciji `Adaptive_PI_Simulation.slx`, pokreće se simulacija sistema upravljanja.
Ova simulacija i rezultati simulacije su opisani u fajlu [Documentation](./Documentation) , gdje su dati i vremenski dijagrami parametara modela, koji su rezultat računanja impementiranog algoritma estimacije. Ovdje je dokumentovano poređenje rezultata za implementirani blok sa blokom *Recursive Polynomial Model Estimator* iz simulinkove biblioteke blokova.



---

## Licenca

Ovaj projekt je licenciran pod MIT licencom. Detalji licence su dati u datoteci [LICENSE](./LICENSE).
