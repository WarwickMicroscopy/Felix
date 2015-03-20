#!/bin/bash

binSIM=felixsim
binREFINE=felixrefine
#binDRAW=felixdraw

inpfile=felix.inp
scafile=felix.sca
ciffile=felix.cif

#Name of Job
#-----------
Job_Name=FelixSim_Test #job name here (FS... or FR...)


#Directory Input
#---------------
#input file directory - point to where all input file(s) are located (at the very least the .sca file)
# NEEDS TO EXIST
inputfiledir=$HOME/project/KTP_felix/Felix/GaA_test #example

# point this to where the felixsim/refine executable is
# NEEDS TO EXIST
binarydir=$HOME/project/KTP_felix/Felix/src

# submission directory - point to what you want the PBS working directory to be 
# Doesn't need to exist
submitdir=/work/e370/e370/ajmh1001/Archer_test_script

# tmp directory to run felix in
# Doesn't need to exist
bashtmpdir=/work/e370/e370/ajmh1001/RUNS/Archer_test_script_temp


#Input Arguments - only input arguments required are the number of nodes and walltime 
#---------------
# settings for parallel submission
nodes=${1:-1}
wtime=${2:-00:01:00} 
JobId=${3:-1}

let cores=${nodes}
let mppwidth=${cores}

#from number of nodes, determine the number of processes
dummyprocesses=$((${cores}*24))
processes=${4:-${dummyprocesses}}

let ranks=${processes}

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


#loop starts here for looped input file parameters or CIFs#
#---------------------------------------------------------#

echo "Submitting jobs for ${cores} cores, ${ranks} MPI ranks, 0 OpenMP threads with maximum wall time ${wtime}"

[ -d ${submitdir} ] || mkdir -p ${submitdir}
cd ${submitdir}
job_dir=${submitdir}

#if looped, change job file name to depend on each iteration
job_file=`printf "${Job_Name}_${ranks}.sh"`

echo "Job file location:" ${job_dir}/${job_file}

[ -d ${job_dir} ] || mkdir -p ${job_dir}


#write PBS scripts
cat > ${job_dir}/${job_file} << EOD
#!/bin/bash --login
##PBS -M a.j.m.hubert@warwick.ac.uk
#PBS -m a
#PBS -r y
#PBS -V
##PBS -k oe
#PBS -j oe

#       The jobname
#PBS -N FelixSim_${ranks}

#       The total number of parallel tasks for your job.
#PBS -l select=${cores}

#       Specify the wall clock time required for your job.
#PBS -l walltime=${wtime}

#       Specify which budget account that your job will be charged to.
#PBS -A e370
  
# Make sure any symbolic links are resolved to absolute path
export PBS_O_WORKDIR=\$(readlink -f \$PBS_O_WORKDIR)

# The base directory is the directory that the job was submitted from.
# All simulations are in subdirectories of this directory.
basedir=\$PBS_O_WORKDIR
echo "basedir=" \${basedir}

# The binary directory is the directory where the scripts are
echo "binarydir (felixsim/ refine location)= " ${binarydir}

# do the run in a TMPdir for safekeeping
tmpdir=${bashtmpdir}/`basename ${job_file} .sh`
[ -d \${tmpdir} ] || mkdir -p \${tmpdir}

# construct the input files

cd \${tmpdir}
echo "Hostname: " $HOSTNAME
echo "Current Working Directory = "
pwd

rm -rf $inpfile

touch $inpfile

echo "# Input file for felixsim version :VERSION: Build :BUILD:"               >> $inpfile
echo "# ------------------------------------"                                  >> $inpfile
echo ""                                                                        >> $inpfile
echo "# ------------------------------------"                                  >> $inpfile
echo "# felixsim input"                                                        >> $inpfile
echo ""                                                                        >> $inpfile
echo "# control flags"                                                         >> $inpfile
echo "IWriteFLAG                = 3"                                           >> $inpfile
echo "IImageFLAG                = 0"                                           >> $inpfile
echo "IOutputFLAG               = 0"                                           >> $inpfile
echo "IBinorTextFLAG            = 0"                                           >> $inpfile 
echo "IScatterFactorMethodFLAG  = 0"                                           >> $inpfile
echo "ICentralBeamFLAG          = 1"                                           >> $inpfile
echo "IMaskFLAG                 = 1"                                           >> $inpfile
echo "IZolzFLAG                 = 1"                                           >> $inpfile
echo "IAbsorbFLAG               = 1"                                           >> $inpfile
echo "IAnisoDebyeWallerFlag     = 0"                                           >> $inpfile
echo "IBeamConvergenceFLAG      = 1"                                           >> $inpfile
echo "IPseudoCubicFLAG          = 0"                                           >> $inpfile
echo "IXDirectionFLAG           = 0"                                           >> $inpfile
echo ""                                                                        >> $inpfile
echo "# radius of the beam in pixels"                                          >> $inpfile
echo "IPixelCount               = 64"                                          >> $inpfile
echo ""                                                                        >> $inpfile
echo "# beam selection criteria"                                               >> $inpfile
echo "IMinReflectionPool        = 2500"                                        >> $inpfile
echo "IMinStrongBeams           = 600"                                         >> $inpfile
echo "IMinWeakBeams             = 0"                                           >> $inpfile
echo "RBSBMax                   = 0.1"                                         >> $inpfile
echo "RBSPMax                   = 0.1"                                         >> $inpfile
echo "RConvergenceTolerance (%) = 1.0"                                         >> $inpfile
echo ""                                                                        >> $inpfile
echo "# crystal settings"                                                      >> $inpfile
echo "RDebyeWallerConstant      = 0.4668"                                      >> $inpfile
echo "RAbsorptionPer            = 2.9"                                         >> $inpfile
echo ""                                                                        >> $inpfile
echo "# microscope settings"                                                   >> $inpfile
echo "ROuterConvergenceAngle    = 10.0"                                        >> $inpfile
echo "RInnerConvergenceAngle    = 0.0"                                         >> $inpfile
echo "IIncidentBeamDirectionX   = 0"                                           >> $inpfile
echo "IIncidentBeamDirectionY   = 1"                                           >> $inpfile
echo "IIncidentBeamDirectionZ   = 0"                                           >> $inpfile
echo "IXDirectionX              = 0"                                           >> $inpfile
echo "IXDirectionY              = -1"                                          >> $inpfile
echo "IXDirectionZ              = 1"                                           >> $inpfile
echo "INormalDirectionX         = 0"                                           >> $inpfile
echo "INormalDirectionY         = 1"                                           >> $inpfile
echo "INormalDirectionZ         = 0"                                           >> $inpfile
echo "RAcceleratingVoltage (kV) = 200.0"                                       >> $inpfile
echo ""                                                                        >> $inpfile
echo "# Image Output Options"                                                  >> $inpfile
echo ""                                                                        >> $inpfile
echo "RInitialThickness        = 1400.0"                                       >> $inpfile
echo "RFinalThickness          = 3000.0"                                       >> $inpfile
echo "RDeltaThickness          = 20.0"                                         >> $inpfile
echo "IReflectOut              = 70"                                           >> $inpfile
echo  ""                                                                       >> $inpfile   
echo "# felixrefine Input"                                                     >> $inpfile
echo ""                                                                        >> $inpfile
echo "#Refinement Specific Flags"                                              >> $inpfile
echo "IRefineModeFLAG = 0"                                                     >> $inpfile
echo "IWeightingFLAG = 0"                                                      >> $inpfile
echo ""                                                                        >> $inpfile
echo "# Debye Waller Factor Iteration"                                         >> $inpfile
echo ""                                                                        >> $inpfile
echo "RInitialDebyeWallerFactor = 0.1"                                         >> $inpfile
echo "RFinalDebyeWallerFactor = 1.0"                                           >> $inpfile
echo "RDeltaDebyeWallerFactor = 0.1"                                           >> $inpfile
echo "IElementsforDWFchange = {0}"                                             >> $inpfile
echo ""                                                                        >> $inpfile
echo "# Ug Iteration"                                                          >> $inpfile
echo ""                                                                        >> $inpfile
echo "INoofUgs = 10"                                                           >> $inpfile
echo "RLowerBoundUgChange = 50.0"                                              >> $inpfile
echo "RUpperBoundUgChange = 50.0"                                              >> $inpfile
echo "RDeltaUgChange = 50.0"                                                   >> $inpfile
echo ""                                                                        >> $inpfile
echo "# Structural Refinement"                                                 >> $inpfile
echo ""                                                                        >> $inpfile
echo "IAtomicsSites = (1,2,3,4,5,6,7)"                                         >> $inpfile
echo ""                                                                        >> $inpfile
echo "# Refinement Output"                                                     >> $inpfile
echo ""                                                                        >> $inpfile
echo "IPrint = 10"                                                             >> $inpfile
echo ""                                                                        >> $inpfile 

cat $inpfile
ls -al \${tmpdir}

#copy input files to temp directory
cp ${inputfiledir}/* \${tmpdir}
echo "copying input files from ${inputfiledir}"
echo "to temp directory: \${tmpdir}"

#if looping over ciffile use this statement to copy each individual one to its tmp folder
#cp \$cifdir/... \${tmpdir}/felix.cif

# Make sure ranks are numbered appropriately
export MPICH_RANK_REORDER_METHOD=1 
export MPICH_PTL_MATCH_OFF=1
export MPICH_FAST_MEMCPY=1
export IPATH_NO_CPUAFFINITY=1

echo "--- starting the SIM run"
aprun -n ${ranks} -N 24 ${binarydir}/felixsim

#echo "---starting the REFINE run"
#aprun -n ${ranks} -N 24 ${binarydir}/felixrefine

echo "--- starting the IMAGE post-processing run"
echo "Output files to be copied"

ls-l

echo "--- starting the image part"

# copy the result of the run in the detination for safekeeping

targetdir=$HOME/RUNS/`basename ${job_file} .sh`
[ -d \${targetdir} ] || mkdir -p \${targetdir}
cp -vr * \${targetdir}

echo "copying output files from" 
pwd
echo "to \${targetdir}" 

wait
#exit 0
EOD

#job submission arguments, uncomment queue required

chmod 755 ${job_dir}/${job_file}
#(cd ${job_dir} ; qsub ./${job_file})
#(cd ${job_dir} ; qsub -q long ./${job_file})
#(cd ${job_dir} ; qsub -q short ./${job_file})
#(cd ${job_dir} ; qsub -q low ./${job_file})
#(cd ${job_dir} ; qsub -q serial ./${job_file})

