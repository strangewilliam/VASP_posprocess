# $La_2CuO_4$

$La: 5d^16s^2$  3

$Cu: 4s^13d^{10}$  11

$O: 2s^22p^4$  6

41

```shell
#!/bin/bash
#La2CuO4_dos_SOC_U_cal

job_name="La2CuO4_dos_SOC_U_cal"
describe="La2CuO4 dos calculation with SOC,with DFT+U, starting from the static calculation, ISYM=2"
start_time=`date`
job_dir=`pwd`
hom1=`pwd`
echo "The home directory is ${hom1}."

######################################
#############  Preparation  #############
######################################
CPU_num=48

methods="PBE"
ISYM_set=2 # consider the symmetry during calculation
SYSTEM_name="La2CuO4"
ENCUT_set=500
ISIF_set=2
ISMEAR_tot_set=-5
ISMEAR_line_set=0
KPOINTS_set="9   9   9"
NBANDS_set=96
ISPIN_set=2
LDAU_set=".TRUE."
## POSCAR_exp
cat >POSCAR<<!
La2 Cu1 O4
1.0
        4.0086321831         0.0000000000         0.0000000000
       -0.0000005338         4.0086321831         0.0000000000
       -2.0043164854         2.0043153360         6.3174697731
   La   Cu    O
    2    1    4
Direct
     0.649622977         0.350376993         0.299246013
     0.350376993         0.649622977         0.700753987
     0.000000000         0.000000000         0.000000000
     0.500000000         0.000000000         0.000000000
     0.000000000         0.500000000         0.000000000
     0.250000000         0.250000000         0.500000000
     0.750000000         0.750000000         0.500000000
!
elements=`sed -n '6,6p' POSCAR`

## POTCAR
rm POTCAR*
for method in ${methods}
do
for ele in ${elements}
do
cat ~/POTCAR/${method}.54/${ele}/POTCAR >> ./POTCAR_${method}
done
## Find the maximum of the ENMAXs as ENCUT
grep 'ENMAX' POTCAR_${method}>temp
row=`awk 'END{print NR}' temp`
m=0
echo ${m}
for r in $( seq ${row} )
do
  t=`awk NR!=${r}'{next}{printf("%.3f",$3)}' temp`
  judge=`echo "${t}>${m}"|bc`
  if [ "${judge}" = "1" ];then
    m=${t}
  fi	
done
rm temp
m=`awk -v m=${m} 'BEGIN{m*=1.5;printf("%.3f",m)}'`
judge=`echo "${m}>${ENCUT_set}"|bc`
if [ "${judge}" = "1" ];then
  ENCUT_set=${m}
fi
done

## KPOINTS_bulk
cat >KPOINTS<<!
Automatic mesh
0
Gamma
${KPOINTS_set}
0.0 0.0 0.0
!

for method in ${methods}
do
mkdir ${method}
cp POTCAR_${method} ${method}/POTCAR
cp POSCAR KPOINTS ./${method}
cd ${method}
######################################
#########  VASP_static_calculation  #########
######################################
# INCAR_stat
cat >INCAR<<!
# Comments
SYSTEM = ${SYSTEM_name}_stat

# I/O
ISTART = 0
ICHARG = 2
LWAVE = .FALSE.
#LCHARG = .FALSE.

# Parallelization: multicore machines linked by a fast network
LPLANE = .TRUE.
NCORE = 6
LSCALU = .FALSE.
NSIM = 4

# Electronic Relaxation
ENCUT = ${ENCUT_set}
NELM = 200
ALGO = Fast
PREC = Accurate
ISMEAR = ${ISMEAR_tot_set}
EDIFF = 1E-7
LREAL = .FALSE. ! Default. More accurate in the reciprocal space than in the real space

ISYM = ${ISYM_set}  ! consider the symmetry

# Polarization
ISPIN = ${ISPIN_set}

# Band
NBANDS = ${NBANDS_set}
!

## SUB_stat
cat >SUB<<!
#BSUB -q e52680v3ib!
#BSUB -n ${CPU_num}
#BSUB -J ${SYSTEM_name}_stat
#BSUB -o out
#BSUB -e err
module load ips/2019u5
mpirun vasp_std
!

if [ "${LDAU_set}" = ".TRUE." ];then
cat >>INCAR<<!

#DFT+U
LDAU      = .TRUE.
LDAUTYPE  = 2            # Default
LDAUL     = 3 2 0        # La d Cu d O p 
LDAUU     = 6 3 0        # specifies the strength of the effective on-site Coulomb interactions.
LDAUJ     = 0.5 0.5 0    # specifies the strength of the effective on-site exchange interactions.
LMAXMIX   = 6            # d-4   f-6
LDAUPRINT = 0            # Default=0, Write occupancy matrix to the OUTCAR file.
!
fi

pwd
date
bsub -K < SUB
grep "reach" OUTCAR | tail -1
tail -14 OUTCAR
echo "The static calulation has finished!"

######################################
#########  VASP_dos_calculation  #########
######################################
mkdir dos
cp POTCAR POSCAR CHGCAR ./dos
cd dos

## KPOINTS_high_symmetry_line
cat >KPOINTS<<!
k-points along high symmetry lines 
192           ! intersections 24*8
Line-mode 
reciprocal 
0.0 0.0 0.0 ! G
1.0 0.0 0.0 ! Z

1.0 0.0 0.0 ! Z
0.5 0.5 0.0 ! X

0.5 0.5 0.0 ! X
0.0 0.0 0.0 ! G

0.0 0.0 0.0 ! G
0.0 0.0 0.5 ! Z

0.0 0.0 0.5 ! Z
0.5 0.0 0.0

0.5 0.0 0.0
0.5 0.0 0.5

0.5 0.0 0.5
0.25 0.25 0.0

0.25 0.25 0.0
0.25 0.25 0.5
!

## INCAR_dos
cat >INCAR<<!
# Comments
SYSTEM = ${SYSTEM_name}_dos

# I/O
ISTART = 0
ICHARG = 11
LWAVE = .FALSE.
LCHARG = .FALSE.

# Parallelization: multicore machines linked by a fast network
LPLANE = .TRUE.
NCORE = 6
LSCALU = .FALSE.
NSIM = 4

# Electronic Relaxation
ENCUT = ${ENCUT_set}
NELM = 200
ALGO = Fast
PREC = Accurate
ISMEAR = ${ISMEAR_line_set}
EDIFF = 1E-7
LREAL = .FALSE. ! Default. More accurate in the reciprocal space than in the real space

ISYM = ${ISYM_set} ! do not consider the symmetry

# Polarization
ISPIN = ${ISPIN_set}

# Band
NBANDS = ${NBANDS_set}
LORBIT = 12
!

## SUB_dos
cat >SUB<<!
#BSUB -q e52680v3ib!
#BSUB -n ${CPU_num}
#BSUB -J ${SYSTEM_name}_dos
#BSUB -o out
#BSUB -e err
module load ips/2019u5
!

case ${ISPIN_set} in
1)
cat >>SUB<<!
mpirun vasp_std
!
;;
2)
cat >>INCAR<<!

#Coupling
LSORBIT    = .TRUE.
GGA_COMPAT = .FALSE.
!
cat >>SUB<<!
mpirun vasp_ncl
!
;;
esac

if [ "${LDAU_set}" = ".TRUE." ];then
cat >>INCAR<<!

#DFT+U
LDAU      = .TRUE.
LDAUTYPE  = 2            # Default
LDAUL     = 3 2 0        # La d Cu d O p 
LDAUU     = 6 3 0        # specifies the strength of the effective on-site Coulomb interactions.
LDAUJ     = 0.5 0.5 0    # specifies the strength of the effective on-site exchange interactions.
LMAXMIX   = 6            # d-4   f-6
LDAUPRINT = 0            # Default=0, Write occupancy matrix to the OUTCAR file.
!
fi

pwd
date
bsub -K < SUB
grep "reach" OUTCAR | tail -1
tail -14 OUTCAR
echo "The dos calculation has finished!" 

cd ${hom1}
done

##############################
########### Record ###########
##############################
end_time=`date` 
cd ~
cat>>script_log<<!
======================================
Job name:      ${job_name}
Job directory: ${job_dir}
Description:   ${describe}
Start time:    ${start_time}
End time:      ${end_time}
If Checked (✓):
======================================
!

exit 0
```





```shell
# !/bin/bash
if [[ -e POTCAR* ]];
then
  rm POTCAR*
fi
if [ ! -e "POSCAR" ];
then
  echo "Please prepare the POSCAR file!"
  exit 1
else
  elements=`sed -n '6,6p' POSCAR`
fi
judge="0"
while [ ${judge} = "0" ]
do
  echo "Please input the type of pseudo-potential: [LDA/PBE]"
  read -n 3 method
  if [ ${method} == "LDA" || ${method} == "PBE" ];
  then
    ${judge} = "1"
  else
    echo "Wrong input!"
    echo "Please input the type of pseudo-potential: [LDA/PBE]"
  fi
done

for ele in ${elements}
do
cat ~/POTCAR/${method}.54/${ele}/POTCAR >> ./POTCAR_${method}
done

## Find the maximum of the ENMAXs
grep 'ENMAX' POTCAR_${method}>temp
row=`awk 'END{print NR}' temp`
m=0
echo ${m}
for r in $( seq ${row} )
do
  t=`awk NR!=${r}'{next}{printf("%.3f",$3)}' temp`
  judge=`echo "${t}>${m}"|bc`
  if [ "${judge}" = "1" ];then
    m=${t}
  fi	
done
rm temp
echo "ENMAX = " ${m}
```

```shell
# !/bin/bash
! Kpoints


exit
```

