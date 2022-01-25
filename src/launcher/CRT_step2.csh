#!/bin/tcsh

if( $#argv != 1 ) then
  echo "USAGE: [script_name].sh [run_num]"
  exit
endif

# step 2: first cut on 20pC on charge, baseline subtraction for each event and creation of TTree with variables of interest
echo .L analysis_CRT.C++"\n"analysis_CRT \*a = new analysis_CRT\( \"run"$1"_CRTNew.root\" \) "\n"a-\>Loop\(\"run"$1"_ana.root\", 0\) | root

echo Output File: run"$1"_ana.root
