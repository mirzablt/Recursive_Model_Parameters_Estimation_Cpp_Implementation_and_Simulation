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

- **ERLSLibrary.cpp**: Funkcije za strukturiranje podataka i dinamičku alokaciju memorije.
- **BFGSAlgorithm.cpp**: Implementacija BFGS algoritma i funkcija linijskog pretraživanja.
- **MatrixLibrary.cpp**: Operacije sa matricama i vektorima.
- **Recursive_Estimation_wrapper.cpp**: Implementacija rekurzivnih algoritama.
- **MPC_controller_parameters.m**: Učitavanje parametara MPC kontrolera.
- **PSAU_BUS.mat**: Definicija dimenzija parametara za višestruke ulaze.
- **Documentation.pdf**: Dokumentacija projekta i tehnički detalji.

---

## Simulacije

### Adaptivno MPC upravljanje

Simulacija koristi ARX model za kontrolu nelinearnog hemijskog reaktora. Implementirani blok rekurzivnog estimatora koristi Kalmanov filter za ažuriranje parametara.

### PI adaptivno upravljanje

PI regulator se podešava na osnovu parametara ARMAX modela. Estimacija parametara se vrši ERLS metodom sa faktorom iščezavanja.

---

## Licenca

Ovaj projekt je licenciran pod **MIT licencom**. Detalji licence nalaze se u datoteci LICENSE.
