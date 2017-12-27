#!/bin/bash 

# Script generates multi-loop samples

###################################################
# Script paths
###################################################

# Path to script directory
SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"
# Root project directory
ROOT_DIR=$SCRIPT_DIR/../..
# Base data directory
DATA_DIR=$ROOT_DIR/ompg
# Path to fragment libraries directory
FRAG_LIB_DIR=$DATA_DIR/frag_libs
# Path to mutations directory
MUT_DIR=$DATA_DIR/mutations
# Path to loop regions directory
LOOP_REGION_DIR=$DATA_DIR/loop_regions
# Path to simulation specifications directory
SIM_INI_DIR=$LOOP_REGION_DIR
# Path to template pdbs
TEMPLATE_DIR=$DATA_DIR/templates
# Base directory containing membrane pdbs
MEMB_DIR=$DATA_DIR/membrane
# Path to base output directory
BASE_OUTPUT_DIR=$DATA_DIR/output/multi_loop
# Path for storing logs
LOG_DIR=$DATA_DIR/logs/multi_loop
# Path for storing energy stats
STAT_DIR=$DATA_DIR/stats/multi_loop
# Path to executable directory
EXE_DIR=$ROOT_DIR
# Name of executable
EXE_NAME=mDisgro

###################################################
# Arrays
###################################################

# Simulation identifiers
SIM_ID=("wt")
# Must be parallel (same size) to SIM_ID
# If i-th value is 0, then SIM_ID[i] does not need mutations applied
# If i-th value is 1, then MUT_DIR/SIM_ID[i].ini identifies the mutations to apply
IS_MUT=(0)

# The template pdb structures to use
# 2iwv -> open conformation crystal structure captured at pH 7
# 2iww -> closed conformation crystal structure captured at pH 5
TEMPLATE=("2iwv" "2iww")

# Whether or not to use loop rotamers
ROT=("-norot" "")

# Overlap factor for adjacent residues (OFA)
OFA=("0.80" "0.85")
# Number of jobs - must be parallel to OFA array
NUM_JOBS=(100 200)

###################################################
# Scalars
###################################################

# Timestamp script execution - append this value to each job prefix
# to avoid overwriting previously generated loops
SCRIPT_ID=$(date +%Y%m%d%H%M%S)

# Number of backbone distance states to sample
NUM_DIST_STATES=32

# Number of side chain distance states to sample
NUM_SC_STATES=16

# Overlap factor for non-adjacent residues
OFNA="0.95"

# Overlap factor for membrane atoms
OFMEMB="0.75"

# Number of attempted samples per job
# => Total number of attempted samples is NUM_JOBS * ATT_PER_JOB
NUM_CONF_PER_JOB="25000"

# Maximum number of backbone clashes allowed
MAX_BB_CLASHES=0

# Maximum number of side chain clashes allowed
MAX_SC_CLASHES=15

# Maximum number of loop restarts when growing a multi-loop sample. A sample
# is tossed and a new sample is attempted from scratch if the restarts
# exceeds this value.
MAX_MU_RESTART=25

# Maximum number of fragments to check when restarting a loop within a
# multi-loop sample
MAX_FRAG_CHECKS=50

# File specifying loop regions
LOOP_REGION_PATH=$LOOP_REGION_DIR/loops_1_2_3_5_6_7.txt

###################################################
# Timestamp utils
###################################################

# From http://stackoverflow.com/questions/17066250/create-timestamp-variable-in-bash-script
timestamp() {
  date +"%T"
}

###################################################
# JOBQUEUE UTILS
###################################################

# From https://pebblesinthesand.wordpress.com/2008/05/22/a-srcipt-for-running-processes-in-parallel-in-bash/

NUM=0
QUEUE=""
MAX_NPROC=12

function queue {
    QUEUE="$QUEUE $1"
    NUM=$(($NUM+1))
}
 
function regeneratequeue {
    OLDREQUEUE=$QUEUE
    QUEUE=""
    NUM=0
    for PID in $OLDREQUEUE
    do
        if [ -d /proc/$PID  ] ; then
            QUEUE="$QUEUE $PID"
            NUM=$(($NUM+1))
        fi
    done
}
 
function checkqueue {
    OLDCHQUEUE=$QUEUE
    for PID in $OLDCHQUEUE
    do
        if [ ! -d /proc/$PID ] ; then
            regeneratequeue # at least one PID has finished
            break
        fi
    done
}

###################################################
# Queue jobs
###################################################

# Iterate over simulation identifiers(i.e. mutations and wild types)
for ix_sim_id in "${!SIM_ID[@]}"
do

    # Name of simulation
    sim_id="${SIM_ID[$ix_sim_id]}"
    # Mutant status: 1 -> mutant, 0 -> wild type
    is_mut="${IS_MUT[$ix_sim_id]}"
    
    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    echo "$(timestamp): simulation: $sim_id"
    echo "$(timestamp): is_mut: $is_mut"
    
    # Iterate over template PDBs
    for template in "${TEMPLATE[@]}"
    do
        # Path to template PDB
        template_pdb_path="$TEMPLATE_DIR/$template.pdb"
        # Path to membrane PDB
        memb_pdb_path="$MEMB_DIR/$template.DMPC.pdb"
        # Prefix for this (simulation, template) pair
        sim_template_prefix="$sim_id"_$template
        # Path to simulation .ini - (ini gives identifiers of frag libs)
        sim_ini_path=$SIM_INI_DIR/sim_"$sim_template_prefix".frag.ini
        # directory to write generated loop pdbs
        pdb_outdir=$BASE_OUTPUT_DIR/$sim_template_prefix/pdb
        
        # Make output directory for this (simulation, template) pair if not existing
        # http://stackoverflow.com/questions/1731767/how-to-create-nonexistent-subdirectories-recursively-using-bash
        #   -p, --parents
        #   no error if existing, make parent directories as needed
        mkdir -p $pdb_outdir
        
        # Directory to write stats for these fragment libraries
        stat_dir="$STAT_DIR/$sim_template_prefix"
        mkdir -p $stat_dir
        
        # Make log directory
        log_dir=$LOG_DIR/$sim_template_prefix/pdb
        mkdir -p $log_dir
        
        echo "+++++++++++++++++++++++++++++"
        echo "$(timestamp): Applying template: $template_pdb_path"
        echo "$(timestamp): Using fragments: $sim_ini_path"
        echo "$(timestamp): Using membrane: $memb_pdb_path"
        echo "$(timestamp): pdb output dir: $pdb_outdir"
        echo "$(timestamp): stats dir: $stat_dir"
        echo "$(timestamp): log dir: $log_dir"
        
         # Iterate over usage of loop rotamers
        for ix_rot in "${!ROT[@]}"
        do
            # Rotamer library usage command line
            rot_cmd="${ROT[$ix_rot]}"
            
            echo "-----------------------------------"
            echo "$(timestamp): rot switch: $rot_cmd"
            
             # Iterate over overlap factor for adjacent loop residues
            for ix_ofa in "${!OFA[@]}"
            do
                # Collision overlap factor for adjacent residues
                ofa="${OFA[$ix_ofa]}"
                # Number of jobs
                num_jobs="${NUM_JOBS[$ix_ofa]}"
                
                echo "***********************************"
                echo "$(timestamp): ofa: $ofa"
                echo "$(timestamp): num jobs: $num_jobs"
                
                # Queue the jobs - C style loop requires double parentheses
                # http://tldp.org/LDP/abs/html/loops1.html
                for ((jid=1; jid <= num_jobs; jid++))
                do
                    # -job_prefix param - uniquely identifies the job that created the .pdb files,
                    # necessary to avoid overwriting .pdbs from other jobs
                    job_prefix="$sim_id"."$template"."rot$ix_rot"."ofa$ofa"."ofna$OFNA"."ofm$OFMEMB"."nds$NUM_DIST_STATES"."nsc$NUM_SC_STATES"."mxbb$MAX_BB_CLASHES"."mxsc$MAX_SC_CLASHES"."$SCRIPT_ID"."$jid"
                    
                    # Path to log file for this job
                    log_path=$log_dir/$job_prefix.txt
                    
                    # Build command line string
                    cmd="./$EXE_NAME"
                    cmd="$cmd -f $template_pdb_path"
                    cmd="$cmd -n $NUM_CONF_PER_JOB -confkeep $NUM_CONF_PER_JOB -pdbout $NUM_CONF_PER_JOB"
                    cmd="$cmd -nds $NUM_DIST_STATES -nscc $NUM_SC_STATES"
                    cmd="$cmd -MaxBBClashes $MAX_BB_CLASHES"
                    cmd="$cmd -MaxSCClashes $MAX_SC_CLASHES"
                    cmd="$cmd -murestart $MAX_MU_RESTART"
                    cmd="$cmd -MaxFragChecks $MAX_FRAG_CHECKS"
                    cmd="$cmd -ofa $ofa -ofna $OFNA -ofmemb $OFMEMB"
                    cmd="$cmd -pdb_outdir $pdb_outdir"
                    cmd="$cmd -stats_outdir $stat_dir"
                    cmd="$cmd -multiloops $LOOP_REGION_PATH"
                    cmd="$cmd -multilooplib_ini $sim_ini_path"
                    cmd="$cmd -multilooplib_dir $FRAG_LIB_DIR"
                    cmd="$cmd -job_prefix $job_prefix"
                    
                    # Append optional switches
                    
                    # Check for mutations
                    if [ "$is_mut" -eq "1" ]; then
                        cmd="$cmd -muts $MUT_DIR/$sim_id.ini"
                    fi # End check for mutations
                    
                    # Check for usage of rotamer library
                    if [ -n "${rot_cmd}" ]; then
                        cmd="$cmd $rot_cmd"
                    fi # End check for usage of rotamer library
                    
                    echo ""
                    echo "$(timestamp): Log path: $log_path"
                    echo "$(timestamp): Running cmd:"
                    echo "$cmd"
                    
                    # Switch to executable dir
                    pushd ./
                    cd $EXE_DIR
                    # Run command in background
                    $cmd > $log_path &
                    # Capture process id of last process ran
                    PID=$!
                    # Store process id
                    queue $PID
                    popd
                    
                    # Delay queues by a few seconds to help ensure unique seeds between jobs
                    sleep 2s
                    
                    # Check if we need to spin until at least one of the jobs finishes
                    while [ $NUM -ge $MAX_NPROC ]; do
                        checkqueue
                        # Sleep for x seconds
                        sleep 10s
                    done # end spin until queue frees up
                    
                done # end iteration over jobs
                
            done # End iteration over adjacent overlap factors
            
        done # End iteration over usage of loop rotamers
        
    done # end iteration over templates
    
done #end iteration over simulation identifiers

###################################################
# Wait on last set of jobs
###################################################

echo "Waiting on final set of jobs to finish..."
wait

echo "Finished."
