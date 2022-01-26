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
