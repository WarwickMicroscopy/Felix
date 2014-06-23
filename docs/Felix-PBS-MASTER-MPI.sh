#!/bin/bash

binSIM=FelixSim
binDRAW=FelixDraw

pytIMAGE=txt2png.py

inpfile=Felix.inp
ciffile=Felix.cif
scafile=Felix.sca

# point this to where FelixSim and FelixDraw are
binarydir=$HOME/D-LACBED/EXE

# point to where Felix.sca is located
scadir=${binarydir}

submitdir=`pwd`
tmpdir=/tmp/

# settings for parallel submission
cores=${1:-1}
ompthreads=${2:-1}
wtime=${3:-00:10:00}

let ranks=${cores}
let nodes=${cores}
let mppwidth=${cores}
curr_dir=`pwd`

if [ ${ranks} -lt 1 ]
then
    ranks=1
fi

if [ ${nodes} -lt 1 ]
then
    nodes=1
fi

if [ ${ranks} -lt 4 ]
then
    taskspernode=${ranks}
else
    taskspernode=4
fi

if [ ${mppwidth} -lt 4 ]
then
    mppwidth=4
fi

if [ ${ompthreads} -lt 1 ]
then
    ompthreads=1
fi

echo "Submitting jobs for ${cores} cores, ${ranks} MPI ranks, ${ompthreads} OpenMP threads with maximum wall time ${wtime}"

[ -d ${submitdir} ] || mkdir ${submitdir}
job_dir=${submitdir}

for nUg in 1:10
do
for pUg in 1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0
do

job_file=`printf "FS_%1i%1i%1i.sh" "${nUg}" "${pUg}"` 

echo ${job_dir}/${job_file}

[ -d ${job_dir} ] || mkdir ${job_dir}

cat > ${job_dir}/${job_file} << EOD
#!/bin/bash --login
#PBS -l pvmem=500mb
##PBS -M yourname@whereyouare.something
#PBS -m a
#PBS -r y
#PBS -V
##PBS -k oe
#PBS -j oe

#       The jobname
#PBS -N Felix_${cores}_${ompthreads}_${ZX}${ZY}${ZZ}

#       The total number of parallel tasks for your job.
#       This is the sum of the number of parallel tasks required by each
#       of the aprun commands you are using. In this example we have
#       ${mppwidth} tasks
##PBS -l mppwidth=${mppwidth}
#PBS -l nodes=${nodes}:ppn=1

#       Specify how many processes per node.
##PBS -l mppnppn=32

#       Specify the wall clock time required for your job.
#PBS -l walltime=${wtime}

#       Specify which budget account that your job will be charged to.
##PBS -A midpluswarwick
  
# Make sure any symbolic links are resolved to absolute path
export PBS_O_WORKDIR=\$(readlink -f \$PBS_O_WORKDIR)

# The base directory is the directory that the job was submitted from.
# All simulations are in subdirectories of this directory.
basedir=\$PBS_O_WORKDIR
echo "basedir=" \${basedir}

# The binary directory is the directory where the scripts are
echo "binarydir=" ${binarydir}

# do the run in a TMPdir for safekeeping
tmpdir=/tmp/`basename ${job_file} .sh`

[ -d \${tmpdir} ] || mkdir \${tmpdir}

cp ${binarydir}/${binSIM} \${tmpdir}/
cp ${binarydir}/${binDRAW} \${tmpdir}/
cp \${basedir}/${ciffile} \${tmpdir}/
cp ${scadir}/${scafile} \${tmpdir}/

# construct the input files

cd \${tmpdir}
echo $HOSTNAME
pwd

touch $inpfile

echo "# Input file for Felix \$Revision: 1.5 $"  > $inpfile
echo "# ------------------------------------"   >> $inpfile
echo ""                                         >> $inpfile
echo "# ------------------------------------"   >> $inpfile
echo "# BLOCH input"                            >> $inpfile
echo ""                                         >> $inpfile
echo "# control flags"                          >> $inpfile
echo "IWriteFLAG                = 2"            >> $inpfile
echo "IImageFLAG                = 1"            >> $inpfile
echo "IOutputFLAG               = 1"            >> $inpfile
echo "IBinorTextFLAG            = 1"            >> $inpfile
echo "IScatterFactorMethodFLAG  = 0"            >> $inpfile
echo "ICentralBeamFLAG          = 1"            >> $inpfile
echo "IMaskFLAG                 = 0"            >> $inpfile
echo "IZOLZFLAG                 = 1"            >> $inpfile
echo "IAbsorbFLAG               = 1"            >> $inpfile
echo "IAnisoDebyeWallerFlag     = 0"            >> $inpfile
echo "IBeamConvergence          = 1"            >> $inpfile
echo "IPseudoCubicFLAG          = 0"            >> $inpfile
echo "IXDirection FLAG          = 0"            >> $inpfile 
echo "IDevFLAG                  = 1"            >> $inpfile     
echo ""                                         >> $inpfile
echo "# radius of the beam in pixels"           >> $inpfile
echo "IPixelCount               = 4"            >> $inpfile
echo ""                                         >> $inpfile
echo "# beam selection criteria"                >> $inpfile
echo "IMinReflectionPool        = 100"          >> $inpfile
echo "IMinStrongBeams           = 7"            >> $inpfile
echo "IMinWeakBeams             = 10"           >> $inpfile
echo "RBSBMax                   = 0.1"          >> $inpfile
echo "RBSPMax                   = 0.1"          >> $inpfile
echo "RConvergenceTolerance (%) = 1"            >> $inpfile
echo ""                                         >> $inpfile
echo "# crystal settings"                       >> $inpfile
echo "RDebyeWallerFactor        = 0.4668"       >> $inpfile
echo "RAbsoprtionPercentage     = 2.9"          >> $inpfile
echo ""                                         >> $inpfile
echo "# microscope settings"                    >> $inpfile
echo "RConvergenceAngle         = 3.0"          >> $inpfile
echo "IIncidentBeamDirectionX   = 1"            >> $inpfile
echo "IIncidentBeamDirectionY   = 1"            >> $inpfile
echo "IIncidentBeamDirectionZ   = 1"            >> $inpfile
echo "IXDirectionX              = 0"            >> $inpfile
echo "IXDirectionY              = 1"            >> $inpfile
echo "IXDirectionZ              = -1"           >> $inpfile
echo "INormalDirectionX         = 1"            >> $inpfile 
echo "INormalDirectionY         = 1"            >> $inpfile
echo "INormalDirectionZ         = 1"            >> $inpfile
echo "RAcceleratingVoltage (kV) = 200.0"        >> $inpfile
echo ""                                         >> $inpfile
echo "# Ug Iteration                            >> $inpfile
echo "INoofUgs                  = ${nUg}"       >> $inpfile
echo "RPercentageUgChange       = ${pUg}"       >> $inpfile
echo ""                                         >> $inpfile
echo "#Refinement Specific Flags"               >> $inpfile
echo "IImageOutputFLAG          = 1"            >> $inpfile
echo ""                                         >> $inpfile
echo "# ------------------------------------"   >> $inpfile
echo "# Draw input"                             >> $inpfile
echo ""                                         >> $inpfile
echo "# sample thickness loop (Angstrom)"       >> $inpfile
echo "RInitialThickness        = 100.0"         >> $inpfile
echo "RFinalThickness          = 1000.0"        >> $inpfile
echo "RDeltaThickness          = 100.0"         >> $inpfile
echo "IReflectOut              = 100"           >> $inpfile

cat $inpfile
ls -al \${tmpdir}

# Make sure ranks are numbered appropriately
export MPICH_RANK_REORDER_METHOD=1 
export MPICH_PTL_MATCH_OFF=1
export MPICH_FAST_MEMCPY=1
export IPATH_NO_CPUAFFINITY=1
  
# Launch the parallel job, ${ranks} MPI ranks in ${taskspernode} tasks per node with ${ompthreads} threads
export OMP_NUM_THREADS=${ompthreads}

echo "--- starting the SIM run"
mpirun -n ${cores} \${tmpdir}/${binSIM} 

echo "--- starting the DRAW run"
mpirun -n ${cores} \${tmpdir}/${binDRAW} 

echo "--- starting the IMAGE post-processing run"
pwd
ls

echo "--- starting the image part"
# copy the result of the run in the detination for safekeeping
targetdir=\${basedir}/`basename ${job_file} .sh`

[ -d \${targetdir} ] || mkdir \${targetdir}

cp -vr * \${targetdir}

wait
#exit 0
EOD

chmod 755 ${job_dir}/${job_file}
#(cd ${job_dir} ; qsub -q serial ./${job_file})
#(cd ${job_dir} ; qsub -q parallel ./${job_file})
#(cd ${job_dir} ; qsub -q devel ./${job_file})
(cd ${job_dir} ; qsub -q taskfarm ./${job_file})

done
done
