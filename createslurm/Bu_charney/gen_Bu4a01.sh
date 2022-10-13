input=("1" "2" "3" "4")
for i in {1..50}
       	#${input[*]}
do 
	cp template.sh Bu4_a01_run$i.sh
	echo "RUNDIR=\$SCRATCH/raytracing3D/arrayjobs/Bu_charney/Bu4_a01/2f_run${i}" >> Bu4_a01_run$i.sh
        echo "mkdir -p \$RUNDIR 
cd \$RUNDIR 
cp \$INPDIR/multraj3D.m . 
cp \$SRCDIR/dispersion3D.m . 
cp \$SRCDIR/sympad_.m . 
cp \$SRCDIR/velocity_ou3D.m . 
cp \$SRCDIR/rhs.m . 
cp \$SRCDIR/savedata.m . 
cp \$SRCDIR/initialize.m . 
cp \$SRCDIR/time_stepping.m . 
echo \"Starting run at: 'date'\"

matlab -nodisplay -r \"\$mfile; exit()\">runlog
echo \"Job finished at: 'date'\"
exit
	" >> Bu4_a01_run$i.sh
	sbatch Bu4_a01_run$i.sh
done

