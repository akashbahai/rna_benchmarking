for i in `find */*remapped.pdb`;do rm $i;done
for i in `find */*temp.pdb`;do rm $i;done
for i in `find */*/*.pdb`;do rm $i;done
for i in `find */*fragments*`;do rm $i;done