  
job_ids=()

#source /opt/intel_2021/oneapi/setvars.sh
#conda activate phonopy

cp ../CONTCAR POSCAR
cp ../INCAR .
sed -i '/ISIF/d; /IBRION/d; /EDIFF/d; /IALGO/d; /LREAL/d; /LWAVE/d; /LCHARGE/d; /PREC/d; /ISMEAR/d; /SIGMA/d;' INCAR

sed -i '$a PREC = Accurate' INCAR
sed -i '$a IBRION = -1' INCAR
sed -i '$a EDIFF = 1e-08' INCAR
sed -i '$a ISMEAR = 0; SIGMA = 0.01' INCAR
sed -i '$a IALGO = 38' INCAR
sed -i '$a LREAL = .FALSE.' INCAR
sed -i '$a LWAVE = .FALSE.' INCAR
sed -i '$a LCHARG = .FALSE.' INCAR

cp ../KPOINTS .
sed -i '$d' KPOINTS
echo "2 2 2" >> KPOINTS

cp ../POTCAR .
cp ../queue_script .

phonopy -d --dim 2 2 2 --pa auto -c POSCAR

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
        job_id=$(qsub queue_script)
        job_ids+=("$job_id")
        echo "Submitted JobID: $job_id" > ~/last_phonopyjob.txt
        cd ..

done

#conda deactivate 
