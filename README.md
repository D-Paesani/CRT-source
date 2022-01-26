# CRT-source

per scaricare la repo e creare tutte le cartelle:


curl https://raw.githubusercontent.com/D-Paesani/CRT-source/main/download_repo.sh  | bash


schema cartelle (da ricrearsi quando si scarica la repo CRT-source):
- CRT-analysis
  - data
    - step1 (symlink alla cartella coi txt)
    - step2
    - template
    - calibration
      - luts_s3
      - luts_s3p
    - step3p
    - step3
  - CRT-source
    - .git
    - launcher
    - includes
    - - questo file README
    - - files .C
  - x



L'analisi (step3) è divisa in 3 parti:
step3 p1 (salva i picchi delle landau e gli offset temporali in data/calibration/luts_s3p e i plot in data/step3p )
step3 p2 (salva gli offset in Z in data/calibration/luts_s3p )
step3 (l'analisi vera e propria)

Per lanciare tutte le 3 parti in sequenza, si può eseguire per esempio:
launcher/CRT_step3all.csh run182

I risultati dello step2 (i vecchi run182_ana.root) vanno in data/step2/run182_s2.root
in data/step1 verrà messo un symlink alla cartella coi txt, semmai dovessimo fare da capo anche lo step2 


cose da fare:
- sistemare il nome histTag, mettendo una funzione virtuale che cambia il formato una tantum per tutta la class, da definire dentro i .C dove si fa il loop. 
Se ci sono histboxes "diversi", si gestisce dentro processObj
- definire due costruttori per TH1 e TH2 (leggere da codice vecchio)
