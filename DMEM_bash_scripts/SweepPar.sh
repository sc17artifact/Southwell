#!/bin/bash

sweep_low=1
sweep_high=10
sweep_incr=1
#delay_arr=($(seq 0 1 10));
delay_arr=2
sweep_params='-sweep_low '$sweep_low' -sweep_high '$sweep_high' -sweep_incr '$sweep_incr

export CSP_NG=1
mat_str=${1}
mat_file="./mat_files/${mat_str}.mtx.bin"
N=${2}
ppn=32
n=$(($ppn*$N))
num_samps=50 #$((($N+3)/4))

data_dir='./data/SweepPar/'${mat_str}'/'
exper_str='SweepPar'
mkdir -p $data_dir
sj_gs_file=$data_dir'sj_gs_'$exper_str'_N'$N'_'$mat_str'.txt'
sps_gs_file=$data_dir'sps_gs_'$exper_str'_N'$N'_'$mat_str'.txt'
sj_direct_file=$data_dir'sj_direct_'$exper_str'_N'$N'_'$mat_str'.txt'
sps_direct_file=$data_dir'sps_direct_'$exper_str'_N'$N'_'$mat_str'.txt'
sds_gs_file=$data_dir'sds_gs_'$exper_str'_N'$N'_'$mat_str'.txt'
sds_direct_file=$data_dir'sds_direct_'$exper_str'_N'$N'_'$mat_str'.txt'
mcgs_file=$data_dir'mcgs_'$exper_str'_N'$N'_'$mat_str'.txt'

#rm -f $sj_gs_file
#rm -f $sj_direct_file
#rm -f $sps_gs_file
rm -f $sps_direct_file
#rm -f $sds_gs_file
rm -f $sds_direct_file
#rm -f $mcgs_file

in_str="srun -N $N -n $n ./DMEM_Southwell -b_file InitGuess/${mat_str}.txt -x_zeros $sweep_params -mat_file $mat_file -format_out -num_samps $num_samps -tol 1e-16"

echo "${mat_str} ${N}"

#echo "SJ, GS"
#$in_str -solver sos_sj  | tee -a $sj_gs_file
#echo "SJ, direct"
#$in_str -solver sos_sj  -loc_solver direct | tee -a $sj_direct_file
#echo "SPS, GS"
#$in_str -solver sos_sps | tee -a $sps_gs_file
echo "SPS, direct"
$in_str -solver sos_sps -loc_solver direct | tee -a $sps_direct_file
#echo "SDS, GS"
#$in_str -solver sos_sds | tee -a $sds_gs_file
echo "SDS, direct"
$in_str -solver sos_sds -loc_solver direct | tee -a $sds_direct_file
#echo "MCGS"
#$in_str -solver sos_mcgs | tee -a $mcgs_file


#for delay in ${delay_arr[@]}; do
#	sds_gs_file=$data_dir'sds_gs_'$exper_str'_N'$N'_'$mat_str'_delay'$delay'.txt'
#        rm -f $sds_gs_file
#        sds_direct_file=$data_dir'sds_direct_'$exper_str'_N'$N'_'$mat_str'_delay'$delay'.txt'
#        rm -f $sds_direct_file
#
#	echo "SDS, GS"
#	$in_str -solver sos_sds | tee -a $sds_gs_file
#	echo "SDS, direct"
#	$in_str -solver sos_sds -loc_solver direct | tee -a $sds_direct_file
#done
