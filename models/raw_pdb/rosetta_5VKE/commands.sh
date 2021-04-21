$ROSETTA_TOOLS/remodel/getBluePrintFromCoords.pl -pdbfile input_files/input.pdb > input_files/new_sequence.remodel
$ROSETTA_BIN/remodel.static.linuxgccrelease  @flag_missing_cterm -overwrite
