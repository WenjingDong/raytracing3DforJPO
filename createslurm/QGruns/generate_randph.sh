for i in {1..1250}
       	#${input[*]}
do 
	cp randph.sh randph_run$i.sh
	echo "RUNDIR=\$SCRATCH/raytracing3D/arrayjobs/QGN100/scaled_rerun/randph_run${i}" >> randph_run$i.sh
        echo "mkdir -p \$RUNDIR 
cd \$RUNDIR 
cp \$INPDIR/multraj3D_scaled.m . 
cp \$INPDIR/psi*.mat . 
cp \$SRCDIR/dispersion3D.m . 
cp \$SRCDIR/sympad_.m . 
cp \$SRCDIR/velocity_ou3D.m . 
echo \"Starting run at: 'date'\"

matlab -nodisplay -r \"\$mfile; exit()\">runlog
echo \"Job finished at: 'date'\"
exit
	" >> randph_run$i.sh
	sbatch randph_run$i.sh
done

