#!/bin/bash 

# Script generate raw (unclustered) back bone fragments

###################################################
# Script paths
###################################################

# Path to script directory
SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"
# Root project directory
ROOT_DIR=$SCRIPT_DIR/../..
# Base data directory
DATA_DIR=$ROOT_DIR/ompg
# Output directories
BASE_PDB_OUTPUT_DIR=$DATA_DIR/output
BASE_FRAG_LIB_PDB_OUTPUT_DIR=$BASE_PDB_OUTPUT_DIR/frag_libs/a_raw_libs
# Stats directories
BASE_STATS_DIR=$DATA_DIR/stats
BASE_FRAG_LIB_STATS_DIR=$BASE_STATS_DIR/frag_libs
# Log directories
BASE_LOG_DIR=$DATA_DIR/logs
BASE_FRAG_LIB_LOG_DIR=$BASE_LOG_DIR/frag_libs/a_raw_libs
# Executable name and directory
EXE_NAME=mDisgro
EXE_DIR=$ROOT_DIR
# Base directory containing template pdbs
TEMPLATE_DIR=$DATA_DIR/templates
# Base directory containing mask files
MASK_DIR=$DATA_DIR/loop_regions
# Base directory containing membrane pdbs
MEMB_DIR=$DATA_DIR/membrane
# Base directory containing mutation files
MUT_DIR=$DATA_DIR/mutations

###################################################
# Arrays
###################################################

# The template pdb structures to use
# 2iwv -> open conformation crystal structure captured at pH 7
# 2iww -> closed conformation crystal structure captured at pH 5
TEMPLATE=("2iwv" "2iww")

# Whether or not to use loop rotamers
ROT=("-norot" "")

# The masks to apply
MASK=("" "loops_1_2_3_5_6_7.txt")
# Overlap factor for adjacent residues (OFA) - must be parallel to MASKS array
# If mask is present, then can use higher OFA
OFA=("0.80" "0.85")
# Number of jobs - must be parallel to MASKS array
# If mask is present, less attempts are needed to generate a target number of samples
NUM_JOBS=(1400 600)

# Start/end region parallel arrays - assumed to match loops_1_2_3_5_6_7.txt
START=("18" "54" "97" "177" "217" "259")
END=("29" "65" "106" "188" "234" "267")

# Family of generated loops
LOOP_ID=("loop_1_wt" "loop_2_wt" "loop_3_wt" "loop_5_wt" "loop_6_wt" "loop_7_wt")
# Mutant status - must be parallel to loop id: 1 -> mutant loop, 0 -> wild type loop
# If mutant loop, the mutation applied is $MUT_DIR/$LOOP_ID.ini
IS_MUT=(0 0 0 0 0 0)
# Index into START/END arrays for each element in LOOP_ID. Must be same size
# as LOOP_ID. Defines which (start,end) regions to use.
REGION=(0 1 2 3 4 5)

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

# Iterate over loop family identifiers(i.e mutations and wild types)
for ix_loop_id in "${!LOOP_ID[@]}"
do
    # Name of loop family
    loop_id="${LOOP_ID[$ix_loop_id]}"
    # Offset for (start, end) intervals
    region="${REGION[$ix_loop_id]}"
    # Start residue identifier
    start="${START[$region]}"
    # End residue identifier
    end="${END[$region]}"
    # Mutant status: 1 -> mutant, 0 -> wild type
    is_mut="${IS_MUT[$ix_loop_id]}"
    
    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    echo "$(timestamp): loop family: $loop_id"
    echo "$(timestamp): start: $start"
    echo "$(timestamp): end: $end"
    echo "$(timestamp): is_mut: $is_mut"
    
    # Iterate over template PDBs
    for template in "${TEMPLATE[@]}"
    do
        # Path to template PDB
        template_pdb_path="$TEMPLATE_DIR/$template.pdb"
        # Path to membrane PDB
        memb_pdb_path="$MEMB_DIR/$template.DMPC.pdb"
        # Directory to write resulting fragment libraries
        pdb_outdir="$BASE_FRAG_LIB_PDB_OUTPUT_DIR/$loop_id/$template"
        mkdir -p $pdb_outdir
        # Directory to write stats for these fragment libraries
        stats_dir="$BASE_FRAG_LIB_STATS_DIR/$loop_id/$template"
        mkdir -p $stats_dir
        # Directory to write logs for these fragment libraries
        log_dir="$BASE_FRAG_LIB_LOG_DIR/$loop_id/$template"
        mkdir -p $log_dir
        
        echo "###################################"
        echo "$(timestamp): template: $template"
        echo "$(timestamp): template pdb: $template_pdb_path"
        echo "$(timestamp): membrane pdb: $memb_pdb_path"
        echo "$(timestamp): pdb output dir: $pdb_outdir"
        echo "$(timestamp): stats dir: $stats_dir"
        echo "$(timestamp): log dir: $log_dir"
        
        # Iterate over usage of loop rotamers
        for ix_rot in "${!ROT[@]}"
        do
            # Rotamer library usage command line
            rot_cmd="${ROT[$ix_rot]}"
            
            echo "-----------------------------------"
            echo "$(timestamp): rot switch: $rot_cmd"
            
            # Iterate over usage of masking other simulated loop regions
            for ix_mask in "${!MASK[@]}"
            do
                # The mask identifier - may be empty which means no mask is used
                mask="${MASK[$ix_mask]}"           
                # Collision overlap factor for adjacent residues
                ofa="${OFA[$ix_mask]}"
                # Number of jobs
                num_jobs="${NUM_JOBS[$ix_mask]}"
                
                echo "***********************************"
                echo "$(timestamp): mask: $mask"
                echo "$(timestamp): ofa: $ofa"
                echo "$(timestamp): num jobs: $num_jobs"
                
                # Queue the jobs - C style loop requires double parentheses
                # http://tldp.org/LDP/abs/html/loops1.html
                for ((jid=1; jid <= num_jobs; jid++))
                do
                    # -job_prefix param - uniquely identifies the job that created the .pdb files,
                    # necessary to avoid overwriting .pdbs from other jobs
                    job_prefix="$loop_id"."$template"."rot$ix_rot"."msk$ix_mask"."ofa$ofa"."ofna$OFNA"."ofm$OFMEMB"."nds$NUM_DIST_STATES"."nsc$NUM_SC_STATES"."mxbb$MAX_BB_CLASHES"."mxsc$MAX_SC_CLASHES"."$SCRIPT_ID"."$jid"
                    
                    # Path to log file for this job
                    log_path=$log_dir/$job_prefix.txt
                    
                    # Build command line string
                    cmd="./$EXE_NAME"
                    cmd="$cmd -f $template_pdb_path"
                    cmd="$cmd -memb $memb_pdb_path"
                    cmd="$cmd -start $start -end $end"
                    cmd="$cmd -n $NUM_CONF_PER_JOB -confkeep $NUM_CONF_PER_JOB -pdbout $NUM_CONF_PER_JOB"
                    cmd="$cmd -nds $NUM_DIST_STATES -nscc $NUM_SC_STATES"
                    cmd="$cmd -MaxBBClashes $MAX_BB_CLASHES"
                    cmd="$cmd -MaxSCClashes $MAX_SC_CLASHES"
                    cmd="$cmd -ofa $ofa -ofna $OFNA -ofmemb $OFMEMB"
                    cmd="$cmd -pdb_outdir $pdb_outdir"
                    cmd="$cmd -stats_outdir $stats_dir"
                    cmd="$cmd -job_prefix $job_prefix"
                    
                    # Append optional switches
                    
                    # Check for mutations
                    if [ "$is_mut" -eq "1" ]; then
                        cmd="$cmd -muts $MUT_DIR/$loop_id.ini"
                    fi # End check for mutations
                    
                    # Check for usage of rotamer library
                    if [ -n "${rot_cmd}" ]; then
                        cmd="$cmd $rot_cmd"
                    fi # End check for usage of rotamer library
                    
                    # Check for enabled mask
                    # http://serverfault.com/questions/7503/how-to-determine-if-a-bash-variable-is-empty
                    if [ -n "${mask}" ]; then
                        cmd="$cmd -mask $MASK_DIR/$mask"
                    fi # End check for enabled mask
                    
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
                        sleep 2s
                    done # end spin until queue frees up
                    
                done # end iteration over jobs
                
            done # End iteration over masks
            
        done # End iteration over usage of loop rotamers
        
    done # End iteration over template PDBs
    
done # End iteration of loop family identifiers

###################################################
# Wait on last set of jobs
###################################################

echo "Waiting on final set of jobs to finish..."
wait

echo "Finished."
