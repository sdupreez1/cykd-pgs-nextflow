#!/usr/bin/env bash
set -euo pipefail

project_name="cykd-pgs-nextflow"

root_dir="/re_gecip/renal/sdupreez/${project_name}"
log_dir="${root_dir}/logs"
old_log_dir="${root_dir}/old_logs"

usage() {
   echo 'USAGE:'
   echo ''
   echo './submit.sh -r|--resume OR -n|--new-run -q <queue_name> [-hkmc]'
   echo ''
   echo '-r, --resume:     resume the previous run of the pipeline'
   echo '-n, --new-run:    begin a new run of the pipeline'
   echo '-q, --queue:      queue to submit pipeline job to (does not affect queues that processes are submitted to)'
   echo ''
   echo '-h, --help:       display this help text'
   echo '-k, --keep-logs:  do not move logs from previous run when resuming'
   echo '-m, --move-main:  only move main logs from previous run (keeps all process logs in logs/ directory)'
   echo '-c, --clean-work: clean all data cached from previous runs'
}

date_time=$(date +%Y-%m-%d_%Hh-%Mm-%Ss)
move-logs() {
   local only_main=false
   local only_main_suffix=""
   if [[ $# > 0 && $1 == --only-main ]]; then
      only_main=true
   fi
   $only_main && local date_time_label=${date_time_label}"_main_only" || local date_time_label=${date_time_label}
   $only_main && local name_search="*_${project_name}.*" || local name_search="*"

   mkdir -p "${old_log_dir}"
   mkdir -p "${log_dir}"
   if [[ -n $(ls -A "${log_dir}") ]]; then
         find "${log_dir}" -mindepth 1 -maxdepth 1 -name "${name_search}" -exec mv {} "${old_log_dir}/${date_time_label}" \;
   fi
}

resume="none"
new_run="none"
keep_all_logs=false
clean_work=false
move_main_logs_only=false
queue="none"
while [[ $# > 0 ]]; do
   case $1 in
      -h|--help)       usage;                     exit 1 ;;
      -r|--resume)     resume=true;               shift 1 ;;
      -n|--new-run)    new_run=true;              shift 1 ;;
      -k|--keep-logs)  keep_all_logs=true;        shift 1 ;;
      -m|--move-main)  move_main_logs_only=true;  shift 1 ;;
      -q|--queue)      queue=$2;                  shift 2 ;;
      -c|--clean-work) clean_work=true;           shift 1 ;;
      -*)              echo "Unknown option: $1"; exit 1 ;;
      *)                                          break ;;
   esac
done

if [[ $resume == "none" && $new_rum == "none" ]]; then
   echo 'ERROR: Must specify one of either -r|--resume or -n|--new-run'
   exit 1
fi

if [[ $queue == "none" ]]; then
   echo 'ERROR: Please provide a queue name with -q <queue_name> or --queue <queue_name> (long queue is advised)'
   exit 1
fi

if [[ $keep_all_logs != true ]];
   $move_main_logs_only && move-logs --only-main || move-logs
fi

command="mkdir -p ${root_dir}/execution_details/${date_time}; module load nextflow/26.04.3-with-plugins"
if [[ $clean_work == true ]]; then
   if [[ $resume == true ]]; then
      echo 'ERROR: Cannot specify both --resume an --clean-work flags'
      exit 1
   else
      command+=" nextflow clean -f;"
   fi
fi
command+=" nextflow run ${project_name}.nf"
$resume && command+=" -resume"

mkdir -p ${root_dir}/misc/pgsc_calc
if [[ ! -e ${root_dir}/misc/pgsc_calc/pgsc_calc_container.sif ]]; then
   echo "Downloading pgsc_calc sif image to ${root_dir}/misc/pgsc_calc..."
   module load singularity/4.1.1
   singularity pull ${root_dir}/misc/pgsc_calc/pgsc_calc_containter.sif docker://docker-remote.artifactory.aws.gel.ac/pgscatalog/pgsc_calc:v2-blob 
fi
if [[ ! -e ${root_dir}/misc/pgsc_calc/pgsc_HGDP+1kGP_v1.tar.tsv ]]; then
   echo "Copying HGDP+1kGP refernce panel to ${root_dir}/misc/pgsc_calc..."
   cp /public_data_resources/pgsc/reference_pannels/pgsc_HGDP+1kGP_v1.tar.tsv ${root_dir}/misc/pgsc_calc 
fi

bsub -q ${queue} \
   -P re_gecip_renal \
   -J ${project_name} \
   -o ${log_dir}/%J_${project_name}.out \
   -e ${log_dir}/%J_${project_name}.err \
   -M 3000 \
   -n 3 \
   -cwd ${root_dir} \
   "${command}"
   
