# POSCAR Selective Dynamics mode 设置

## posmove

```shell
# !/bin/bash
# Need POSCAR file and move.in file

#############################################
################  Welcome ###################
#############################################

echo
echo "Welcome!"
echo "This script is used to customize which atoms provided in the POSCAR file are allowed to move."
echo "  You need to provide the POSCAR and the move.in file to generate the POSCAR file for partial structural optimization. The output POSCAR will be named as POSCAR_move."
echo "An example of the move.in file:"
echo "======================================"
printf "%3d%4s\n" "2" "xy"
printf "%3d%4s\n" "4" "xyz"
printf "%3d%4s\n" "5" "z"
printf "%3d%4s\n" "3" "x"
printf "%3d%4s\n" "18" "y"
printf "%3d%4s\n" "30" "z"
echo "======================================"
echo "  On each line, the number i on the left corresponds to the i-th atom provided in the POSCAR file, and the string on the right contains the directions allowed to move."
echo "P.S. If there are repeated numbers, the last one is effective."
read -p "Continue?[y/n]" a
if [ "${a}" != "y" ];then
  echo "Quit."
  echo "Finished."
  echo "Have a good day!"
  exit 0
fi
#############################################
################ Preparation #################
##############################################
# Judge if the POSCAR file exists.
judge=0
if [ ! -e "POSCAR" ];then
  echo "Please provide the POSCAR file!"
  exit 0
fi

# Judge if the POSCAR file contains "Selective dynamics"
selective=0
j=`awk '{if (NR==8){a=toupper(substr($1,1,1))}}END{print a}' POSCAR`
if  [ "${j}" == "S" ];then
  selective=1
fi

# Judge if the POSCAR file is correct.
num_true=`awk -v s=${selective} 'END{print NR-8-s}' POSCAR` 
num_system=`awk 'BEGIN{a=0}{if(NR==7){a=$1+$2+$3+$4}}END{print a}' POSCAR`
if [ "${num_true}" != "${num_system}" ];then
  let judge++
  echo "${judge}. The number of the atoms provided in the POSCAR file is wrong! Please provide the correct POSCAR file!"
fi

# Judge if the move.in file exists
if  [ ! -e "move.in" ];then
  let judge++
  echo "${judge}. Please provide the move.in file!"
fi

if [ "${judge}" != "0" ];then
    echo "Have a good day!"
    exit 0
fi

#############################################
#################  Input  ###################
##############################################

num_atom=`awk 'END{print NR}' move.in`
echo "The number of the atoms in the system: ${num_system}"
echo "The number of atoms to be customized: ${num_atom}"
declare -a num direc num_direc move
num=();direc=();num_direc=();move=();
printf "%12s%17s\n" "Number" "Direction"
for i in $(seq 1 ${num_system})
do
  move[${i}]=" F F F";
done
for i in $(seq 1 ${num_atom})
do
  temp=`awk -v n=${i} '{if (NR==n) {print $1}}' move.in`
  num[${i}]=${temp}
  if (( ${num[${i}]} < 1 )) || (( ${num[${i}]} > ${num_system} ));then
    let judge++
    echo "${judge}. The index of atom in move.in file is wrong!"
    echo "Have a good day!"
    exit 0
  fi
  temp=`awk -v n=${i} '{if (NR==n) {print $2}}' move.in`
  direc[${i}]=${temp}
  num_direc=`echo ${temp}|awk -F"[x-z]" 'END{print NF-1}'`
  x_axis=0;y_axix=0;z_axis=0;
  x_axis=`echo ${temp}|awk -F"x" 'END{print NF-1}'`;
  y_axis=`echo ${temp}|awk -F"y" 'END{print NF-1}'`;
  z_axis=`echo ${temp}|awk -F"z" 'END{print NF-1}'`;
  let "num_direc_temp=x_axis+y_axis*2+z_axis*4"
  case ${num_direc_temp} in
    1) # x=1
      move[${num[${i}]}]=" T F F" 
    ;;
    2) # y=2
      move[${num[${i}]}]=" F T F" 
    ;;
    3) # xy=1+2=3
      move[${num[${i}]}]=" T T F" 
    ;;
    4) # z=4
      move[${num[${i}]}]=" F F T"
    ;;
    5) # xz=1+4=5
      move[${num[${i}]}]=" T F T" 
    ;;
    6) # yz=2+4=6
      move[${num[${i}]}]=" F T T"
    ;;
    7) # xyz=1+2+4=7
      move[${num[${i}]}]=" T T T"
     ;;
    *) # The move.in file is wrong
      let judge++
      echo "${judge}. The move.in file is wrong!"
      echo "Have a good day!"
      exit 0
    ;;
  esac
  num_direc[${i}]=${num_direc_temp}
  printf "%12d%17s\n" ${num[${i}]} ${direc[${i}]}
done

#############################################
#################  Output  ###################
##############################################
# 20.16
case ${selective} in
  0)
    head -8  POSCAR > POSCAR_move
  ;;
  1)
    head -7 POSCAR > POSCAR_move
    head -9 POSCAR|tail -1 >> POSCAR_move
  ;;
  *)
    let judge++
    echo "${judge}. An internal mistake!"
    echo "Have a good day!"
    exit 0
  ;;
esac
n=`awk 'END{print NR}' POSCAR`
j=0
let start=9+selective
for i in $(seq ${start} ${n})
do
  temp=`awk -v j=${i} '{if(NR==j){printf ("%23.16f%23.16f%23.16f", $1,$2,$3)}}' POSCAR`
  let j=i-8-selective
  printf "%s69%7s\n" "${temp}" "${move[${j}]}" >> POSCAR_move
done
sed -i '8iSelective dynamics' POSCAR_move
echo "Congratulations! The POSCAR_move file has been generated successfully ."
echo -e "Do you want to rename the POSCAR_move file as POSCAR?\nType [y] to rename the POSCAR file or any other key to keep the name or [Ctrl+c] to quit."
read judge
echo
if [ "${judge}" = "y" ];then
  mv POSCAR_move POSCAR
fi
echo "POSCAR file"
echo "==========================================="
cat POSCAR
echo "==========================================="
echo
echo "Finished."
echo "Have a good day!"
exit 0
```

## posfix

```shell
# !/bin/bash
# Need POSCAR file and fix.in file

#############################################
################  Welcome ###################
#############################################

echo
echo "Welcome!"
echo "This script is used to customize which atoms provided in the POSCAR file are allowed to fix."
echo "  You need to provide the POSCAR and the fix.in file to generate the POSCAR file for partial structural optimization. The output POSCAR will be named as POSCAR_fix."
echo "An example of the fix.in file:"
echo "======================================"
printf "%3d%4s\n" "2" "xy"
printf "%3d%4s\n" "4" "xyz"
printf "%3d%4s\n" "5" "z"
printf "%3d%4s\n" "3" "x"
printf "%3d%4s\n" "18" "y"
printf "%3d%4s\n" "30" "z"
echo "======================================"
echo "  On each line, the number i on the left corresponds to the i-th atom provided in the POSCAR file, and the string on the right contains the directions allowed to fix."
echo "P.S. If there are repeated numbers, the last one is effective."
read -p "Continue?[y/n]" a
if [ "${a}" != "y" ];then
  echo "Quit."
  echo "Finished."
  echo "Have a good day!"
  exit 0
fi
#############################################
################ Preparation #################
##############################################
# Judge if the POSCAR file exists.
judge=0
if [ ! -e "POSCAR" ];then
  echo "Please provide the POSCAR file!"
  exit 0
fi

# Judge if the POSCAR file contains "Selective dynamics"
selective=0
j=`awk '{if (NR==8){a=toupper(substr($1,1,1))}}END{print a}' POSCAR`
if  [ "${j}" == "S" ];then
  selective=1
fi

# Judge if the POSCAR file is correct.
num_true=`awk -v s=${selective} 'END{print NR-8-s}' POSCAR` 
num_system=`awk 'BEGIN{a=0}{if(NR==7){a=$1+$2+$3+$4}}END{print a}' POSCAR`
if [ "${num_true}" != "${num_system}" ];then
  let judge++
  echo "${judge}. The number of the atoms provided in the POSCAR file is wrong! Please provide the correct POSCAR file!"
fi

# Judge if the fix.in file exists
if  [ ! -e "fix.in" ];then
  let judge++
  echo "${judge}. Please provide the fix.in file!"
fi

if [ "${judge}" != "0" ];then
    echo "Have a good day!"
    exit 0
fi

#############################################
#################  Input  ###################
##############################################

num_atom=`awk 'END{print NR}' fix.in`
echo "The number of the atoms in the system: ${num_system}"
echo "The number of atoms to be customized: ${num_atom}"
declare -a num direc num_direc fix
num=();direc=();num_direc=();fix=();
printf "%12s%17s\n" "Number" "Direction"
for i in $(seq 1 ${num_system})
do
  fix[${i}]=" T T T";
done
for i in $(seq 1 ${num_atom})
do
  temp=`awk -v n=${i} '{if (NR==n) {print $1}}' fix.in`
  num[${i}]=${temp}
  if (( ${num[${i}]} < 1 )) || (( ${num[${i}]} > ${num_system} ));then
    let judge++
    echo "${judge}. The index of atom in fix.in file is wrong!"
    echo "Have a good day!"
    exit 0
  fi
  temp=`awk -v n=${i} '{if (NR==n) {print $2}}' fix.in`
  direc[${i}]=${temp}
  num_direc=`echo ${temp}|awk -F"[x-z]" 'END{print NF-1}'`
  x_axis=0;y_axix=0;z_axis=0;
  x_axis=`echo ${temp}|awk -F"x" 'END{print NF-1}'`;
  y_axis=`echo ${temp}|awk -F"y" 'END{print NF-1}'`;
  z_axis=`echo ${temp}|awk -F"z" 'END{print NF-1}'`;
  let "num_direc_temp=x_axis+y_axis*2+z_axis*4"
  case ${num_direc_temp} in
    1) # x=1
      fix[${num[${i}]}]=" F T sT" 
    ;;
    2) # y=2
      fix[${num[${i}]}]=" T F T" 
    ;;
    3) # xy=1+2=3
      fix[${num[${i}]}]=" F F T" 
    ;;
    4) # z=4
      fix[${num[${i}]}]=" T T F"
    ;;
    5) # xz=1+4=5
      fix[${num[${i}]}]=" F T F" 
    ;;
    6) # yz=2+4=6
      fix[${num[${i}]}]=" T F F"
    ;;
    7) # xyz=1+2+4=7
      fix[${num[${i}]}]=" F F F"
     ;;
    *) # The fix.in file is wrong
      let judge++
      echo "${judge}. The fix.in file is wrong!"
      echo "Have a good day!"
      exit 0
    ;;
  esac
  num_direc[${i}]=${num_direc_temp}
  printf "%12d%17s\n" ${num[${i}]} ${direc[${i}]}
done

#############################################
#################  Output  ###################
##############################################
# 20.16
case ${selective} in
  0)
    head -8  POSCAR > POSCAR_fix
  ;;
  1)
    head -7 POSCAR > POSCAR_fix
    head -9 POSCAR|tail -1 >> POSCAR_fix
  ;;
  *)
    let judge++
    echo "${judge}. An internal mistake!"
    echo "Have a good day!"
    exit 0
  ;;
esac
n=`awk 'END{print NR}' POSCAR`
j=0
let start=9+selective
for i in $(seq ${start} ${n})
do
  temp=`awk -v j=${i} '{if(NR==j){printf ("%23.16f%23.16f%23.16f", $1,$2,$3)}}' POSCAR`
  let j=i-8-selective
  printf "%s69%7s\n" "${temp}" "${fix[${j}]}" >> POSCAR_fix
done
sed -i '8iSelective dynamics' POSCAR_fix
echo "Congratulations! The POSCAR_fix file has been generated successfully ."
echo "POSCAR file"
echo "==========================================="
cat POSCAR
echo "==========================================="
echo
echo "Finished."
echo "Have a good day!"
exit 0
```

