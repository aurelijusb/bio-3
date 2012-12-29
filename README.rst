Bioinfomratics
==============

Compiling/running
-----------------

    python HpvAnalyser.py
    
Libraries and external tools
----------------------------

 * BioPython http://biopython.org/wiki/Main_Page
 * CD-HIT http://weizhong-lab.ucsd.edu/cd-hit/
 * MAfft http://mafft.cbrc.jp/alignment/software/

TODO
----

Užduotis

 Sukurti pobes parinkimo diagnostinei sistemai programą,
skirtą sukurti probe tinkama tik tam tikrų virusų DNR
detekcijai, kuomet padauginimui naudojma PGR (Polimerazinę
grandininę reakciją)
 Sukurtą programą pritaikyti
projektuojant diagnostinę sistemą žmogaus papilomos viruso
(HPV) 16,18,31,33,35,51,52 tipų (sukelia gimdos kaklelio vėžį)
infekcijai nustatyti.
 Diagnostinė sistema turi būti specifiška
ir turi nebūti signalo nepavojingų HPV infekcijų atvejais
(tipai 6, 11, 40, 42, 43, 44, 57, 81).
 Dokumentuota programa turėtų parinkinėti DNR probe L1 geno regione.
 Atitinkama šio geno seka iš HPV 16 tipo genomo pateikta šio
dokumento pabaigoje.
 Šią seką galite naudoti per blast užklausas rinkdami sekas analizei.
 Kuriant algoritmą tikslas būtų sukurti jį tokį, kad būtų parenkamas
kuo mažesnis probiu kiekis tinkantis visų žinomų didelės rizikos
HPV tipų variantų diagnostikai.

II dalies uzduotis:
Užduotis

Parinkti probių sistemą, kurioje būtų kuo mažiau probių, taip kad:
1. Visų probių sekų lydymosi temperatūra turėtų nesiskirti daugiau nei 5 laipsniais (rekomenduojamas ilgis 15-30 bp)
2. Kiekvienai didelio pavojingumo papilomos viruso sekai turėtų būti bent po vieną probę jų rinkinyje, kurios lydymosi temperatūra ant viruso sekos būtų ne mažesnė nei 60 laipsnių (virusas bus detektuojamas).
3. Visų probių lydymosi temperatūra ant   mažo pavojingumo viruso sekos neturėtų būti daugiau nei 55 laipsniai.
4. Regionas, iš kurio parenkamos probės neturėtų būti ilgesnis nei 60bp. Rekomenduojamas probių ilgis 15-30 bp (virusas bus nedetektuojamas).

Lydymosi temperatūrai apskaičiuoti rekomenduoju programą http://sourceforge.net/projects/melting/.
Pradėti rekomenduoju nuo funkcijos, kurios parametrai būtų dvi sekos: prilydomą seka ir seką ant kurios bus lydoma ir ji paskaičiuotų lydymosi temperatūrą paleisdama melting programą.
O poto jūsų fantazijos reikalas.

Galbūt daugeliui lydymosi temperatūra per daug komplikuoja suvokimą.
Pateiku užduotį be jos.

Parinkti probių sistemą, kurioje būtų kuo mažiau probių, taip kad:
1. Visų probių ilgis būtų 30±5 bazių porų (bp).
2. Kiekvienai didelio pavojingumo papilomos viruso sekai turėtų būti bent po vieną probę rinkinyje, taip kad tarp jos ir viruso sekos būtų NEdaugiau nei 2 nesutapimai.
3. Visos probės turėtų bent po 3 nesutapimus su VISOMIS nepavojingų virusų sekomis.
4. Regionas, iš kurio parenkamos probės neturėtų būti ilgesnis nei 60 bp. 



References
----------

 * http://biopython.org/DIST/docs/tutorial/Tutorial.html
 * http://bips.u-strasbg.fr/fr/Tutorials/Comparison/Blast/blastall.html
 * http://ghr.nlm.nih.gov/glossary=probe
 * http://www.ncbi.nlm.nih.gov/nuccore/
