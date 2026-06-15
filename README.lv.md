# 2D Granulārās Gāzes Simulācija

Interaktīva divdimensiju granulārās gāzes simulācija ar neelastīgām sadursmēm, reāllaika vizualizāciju un kinētiskās enerģijas analīzi.

## Pārskats

Šajā projektā ir izstrādāts skaitlisks granulārās gāzes modelis, kurā identiskas daļiņas pārvietojas slēgtā divdimensiju telpā un savstarpēji saduras. Atšķirībā no ideālas gāzes, daļiņu sadursmes ir neelastīgas, tāpēc katras sadursmes laikā daļa kinētiskās enerģijas tiek disipēta.

Simulācija ļauj pētīt:

* kinētiskās enerģijas disipāciju;
* restitūcijas koeficienta ietekmi uz sistēmas dinamiku;
* daļiņu klasterizāciju un neelastīgo kolapsu;
* telpiskās sadalīšanas algoritmu veiktspēju.

Projekts izstrādāts kā vidusskolas pētnieciskais darbs, apvienojot klasisko mehāniku, skaitlisko modelēšanu un zinātnisko programmēšanu Python vidē.

---

## Galvenās iespējas

* 2D cieto sfēru (Hard Sphere) modelis;
* regulējams restitūcijas koeficients (`0.50 ≤ e ≤ 1.00`);
* daļiņu kustības animācija reāllaikā;
* normalizētas kinētiskās enerģijas grafiks;
* ārējās enerģijas pievadīšana ar funkciju **Sakratīt**;
* vairāku eksperimentu rezultātu saglabāšana un salīdzināšana;
* sadursmju meklēšanas optimizācija, izmantojot **Cell List (Spatial Hashing)** algoritmu;
* iebūvēta fizikas modeļa pašpārbaude.

---

## Fizikālais modelis

Simulācijā tiek pieņemti šādi pieņēmumi:

* visas daļiņas ir identiskas un tām ir vienāda masa un rādiuss;
* gravitācijas ietekme netiek ņemta vērā;
* kastes sienas ir absolūti elastīgas;
* daļiņas mijiedarbojas tikai tieša kontakta brīdī;
* neelastīgās sadursmes raksturo restitūcijas koeficients `e`.

Ja `e = 1.0`, sistēmas kopējā kinētiskā enerģija saglabājas konstanta.

Ja `e < 1.0`, sistēma pakāpeniski zaudē kinētisko enerģiju, un noteiktos apstākļos var veidoties daļiņu klasteri un iestāties neelastīgais kolapss.

---

## Izmantotie algoritmi

### Sadursmju detektēšana

Naivā pieeja pieprasa pārbaudīt attālumu starp visiem iespējamiem daļiņu pāriem, kā rezultātā algoritma laika sarežģītība ir:

**O(N²)**

Lai paātrinātu aprēķinus, simulācijas telpa tiek sadalīta regulārā režģī (Cell List metode). Katra daļiņa pārbauda sadursmes tikai savā un blakus esošajās šūnās, tādējādi vidējā gadījumā panākot gandrīz lineāru sarežģītību:

**O(N)**

pie ierobežota daļiņu blīvuma.

---

## Izmantotās tehnoloģijas

* Python 3
* NumPy
* Matplotlib

---

## Instalēšana

Repozitorija klonēšana:

```bash
git clone https://github.com/lietotajvards/granular-gas-simulation.git
cd granular-gas-simulation
```

Atkarību instalēšana:

```bash
pip install numpy matplotlib
```

---

## Programmas palaišana

Grafiskās lietotnes palaišana:

```bash
python granular_gas_simulation.py
```

Fizikas kodola pašpārbaude:

```bash
python granular_gas_simulation.py --self-test
```

Palaišana ar pielāgotiem parametriem:

```bash
python granular_gas_simulation.py --particles 500 --radius 0.007 --restitution 0.9
```

---

## Vadības elementi

| Elements             | Funkcija                                    |
| -------------------- | ------------------------------------------- |
| Restitūcijas slīdnis | Maina restitūcijas koeficientu reāllaikā    |
| Soļi kadrā           | Maina simulācijas ātrumu                    |
| Sakratīt             | Pievada sistēmai papildu kinētisko enerģiju |
| Saglabāt līkni       | Saglabā pašreizējo enerģijas grafiku        |
| Notīrīt              | Dzēš saglabātās līknes                      |

---

## Pielietojums

Simulāciju iespējams izmantot:

* granulāro materiālu dinamikas pētījumiem;
* disipējošu sistēmu modelēšanai;
* statistiskās fizikas un termodinamikas apguvei;
* skaitļošanas fizikas un algoritmu optimizācijas demonstrēšanai;
* zinātniskās vizualizācijas mērķiem.

---

## Autors

Izstrādāts kā vidusskolas pētnieciskais darbs par granulāro gāzu skaitlisko modelēšanu un kinētiskās enerģijas disipācijas procesiem neelastīgu daļiņu sistēmās.
