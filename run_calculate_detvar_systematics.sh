#!/bin/bash

./calculate_detvar_systematics_final.exe -v LYDown -o systematics/detvar/calculate_detvar_systematics_LYDown_ccnueChi2BinningNoBkgHack_ccnumuNewOverflowBinning_output.root -m -c
./calculate_detvar_systematics_final.exe -v LYRayleigh -o systematics/detvar/calculate_detvar_systematics_LYRayleigh_ccnueChi2BinningNoBkgHack_ccnumuNewOverflowBinning_output.root -m -c
./calculate_detvar_systematics_final.exe -v LYAtt -o systematics/detvar/calculate_detvar_systematics_LYAtt_ccnueChi2BinningNoBkgHack_ccnumuNewOverflowBinning_output.root -m -c
./calculate_detvar_systematics_final.exe -v SCE -o systematics/detvar/calculate_detvar_systematics_SCE_ccnueChi2BinningNoBkgHack_ccnumuNewOverflowBinning_output.root -m -c
./calculate_detvar_systematics_final.exe -v recomb2 -o systematics/detvar/calculate_detvar_systematics_recomb2_ccnueChi2BinningNoBkgHack_ccnumuNewOverflowBinning_output.root -m -c
./calculate_detvar_systematics_final.exe -v wiremodX -o systematics/detvar/calculate_detvar_systematics_wiremodX_ccnueChi2BinningNoBkgHack_ccnumuNewOverflowBinning_output.root -m -c
./calculate_detvar_systematics_final.exe -v wiremodYZ -o systematics/detvar/calculate_detvar_systematics_wiremodYZ_ccnueChi2BinningNoBkgHack_ccnumuNewOverflowBinning_output.root -m -c
./calculate_detvar_systematics_final.exe -v wiremodThetaXZ -o systematics/detvar/calculate_detvar_systematics_wiremodThetaXZ_ccnueChi2BinningNoBkgHack_ccnumuNewOverflowBinning_output.root -m -c
./calculate_detvar_systematics_final.exe -v wiremodThetaYZ -o systematics/detvar/calculate_detvar_systematics_wiremodThetaYZ_ccnueChi2BinningNoBkgHack_ccnumuNewOverflowBinning_output.root -m -c