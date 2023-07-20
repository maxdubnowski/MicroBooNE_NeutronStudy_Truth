# MicroBooNE_NeutronStudy_Truth




## Initialize system

source SetUp_FlatTree.sh


## Run through events, apply selection cuts, and fill histograms

root -l script_MicroBooNE_Truth.root


## Change selection cuts, or add plots

emacs NeutrinoSelection.C


## Plot 1D Histograms

root -l PlotRoot.cpp


## Plot 2D Histograms

root -l PlotRoot2D.cpp
