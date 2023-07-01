#phonopy -d --dim 2 2 2 --pa auto -c POSCAR
  
for m in `ls POSCAR-0*`
do
fol=`echo $m|cut -d "-" -f 2`
mkdir disp_$fol
cp POSCAR-$fol disp_$fol/POSCAR
cp INCAR disp_$fol/INCAR
cp POTCAR disp_$fol/POTCAR
cp KPOINTS disp_$fol/KPOINTS
cp queue_script disp_$fol/queue_script

cd disp_$fol/
qsub queue_script
cd ..

done
