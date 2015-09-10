#PBS -l nodes=18:ppn=8
#PBS -l walltime=336:00:00
#PBS -j oe

# change the current working directory to the directory where
# the executable file 'hello' can be found
cd $PBS_O_WORKDIR
echo $PBS_O_WORKDIR

# run the executable file 'hello' using the qmpirun script
mpirun ./../src/bin/psrjmh_ais_auto > ./psrjmh_ais_auto.log
