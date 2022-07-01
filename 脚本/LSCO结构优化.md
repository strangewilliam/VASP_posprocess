

**全程计算，POTCAR选用LDA**

## 结构优化

### 第一步

SUB

```shell
#BSUB -q e52680v3ib!
#BSUB -n 24
#BSUB -J 1_opt
#BSUB -o out
#BSUB -e err
module load ips/2019u5
mpirun vasp_std > log
```

KPOINTS

````shell
Automatic mesh
0
Gamma
2   2   2
0.0 0.0 0.0
````

INCAR

```shell
# Comments
SYSTEM = opt_1

# I/O
ISTART = 0
#LWAVE = .FALSE.
LCHARG = .FALSE.

# Electronic Relaxation
ENCUT = 600
NELM = 100
PREC = Normal
ISMEAR = 0
SIGMA = 0.1
EDIFF = 1E-6
LMAXMIX = 4
LREAL = Auto
ISYM = 0

# Ionic Relaxation
IBRION = 2
ISIF = 2
EDIFFG = 1E-3
NSW = 200

# Bands
NBANDS = 192
#LORBIT = 12

# Parallelization
NPAR = 1
```

https://baijiahao.baidu.com/s?id=1636837921100531467&wfr=spider&for=pc

记录一下计算的时间以及熵值，以及是否收敛（reach关键字）

EENTRO：OUTCAR 中最后一个 EENTRO 值。为最后给出的体系的 “熵” 值。

半导体：ISMEAR=-5

金属：ISMEAR=1，SIGAM=0.1。

可以继续使用 ISMEAR=0，但是需要根据 EENTRO 设置合理的 SIGMA 值。

对于金属，根据 step1 中最后的 EENTRO 值调整 SIGMA 值，尽可能使得 EENTRO/atom 在 1meV ~ 2meV。(SIGMA 越大，EENTRO 最大)。**按照手册，尽量1 meV/atom以内。**

大致判断一下体系是半导体还是金属，进而修改下一步ISMEAR的值，配合着修改SIGMA的值。

https://www.sohu.com/a/320734304_120170422

三. 如何快速判断体系是半导体还是金属

我个人通常习惯做分步结构优化。在第一步结构优化中，使用较低精度，同时设置 ISMEAR = 0 + SIGMA = 0.1，查看结构优化最后在 OUTCAR 中给出的 EENTRO 值。然后通过 EENTRO 值除以体系原子数去判断体系是半导体和金属：

\1) ( EENTRO / 原子数 ) 大于 1meV，体系为金属。

\2) ( EENTRO / 原子数 ) 小于 0.1meV，体系为半导体。这个值非常小时，比如 0.000001，体系带隙一般较大。

\3) ( EENTRO / 原子数 ) 介于 1meV 和 0.1meV，目前还没有具体测试过，不确定。

这一判断的依据是：对于金属体系，ISMAER = 0 会使用 Gaussian smearing，导致在费米面处出现分数占据，进而导致 EENTRO 值不为 0 。而对于较大带隙的半导体，SIGMA = 0.1 还无法导致费米面处出现分数占据，EENTRO 值保持为 0（若增大 SIGMA 到一个很大的值时，EENTRO 也将不再保持为 0，但此时的计算会存在问题）。

这一判断方法是我个人的经验，目前在使用过程中大部分情况下都能给出正确的结果，但不保证对所有体系所有情况都可用。比如在弛豫过程中晶格变化非常大时，用于自洽计算对应的平面波截断球已经严重变形，可能会导致离子步的自洽收敛到的状态本身就有很大的误差。此时再在有误差的自洽结果上用这一方法进行判断是否可行，我不确定。不过解决办法也简单，再用弛豫的结构做一次弛豫就可以使用这一方法进行判断了。

结构优化脚本opt1

```shell
#! /bin/bash

cat >INCAR<<!
# Comments
SYSTEM = opt_1

# I/O
ISTART = 0
#LWAVE = .FALSE.
LCHARG = .FALSE.

# Electronic Relaxation
ENCUT = 600
NELM = 100
PREC = Normal
ISMEAR = 0
SIGMA = 0.1
EDIFF = 1E-6
LMAXMIX = 4
LREAL = Auto
ISYM = 0

# Ionic Relaxation
IBRION = 2
ISIF = 2
EDIFFG = 1E-3
NSW = 200

# Bands
NBANDS = 192
#LORBIT = 12

# Parallelization
NPAR = 1
!

cat >KPOINTS<<!
Automatic mesh
0
Gamma
2   2   2
0.0 0.0 0.0
!

cat >SUB<<!
#BSUB -q e52680v3ib!
#BSUB -n 24
#BSUB -J 1_opt
#BSUB -o out
#BSUB -e err
module load ips/2019u5
mpirun vasp_std > log
!

if [ -e "POSCAR" ];then
  potmake
  typeset -l judge
  read -p "Restrict some atoms[y?]" judge
  if [ "${judge}" == "y" ];then
    echo
    read -p "Allowed to move or fix?[m/f]" judge
    case ${judge} in
      "m")
        posmove
      ;;
      "f")
        posfix
      ;;
      *)
        echo "No restriction added to POSCAR file."
      ;;
     esac
  else
    echo
    echo "No restriction added to POSCAR file."
   fi 
fi

echo
echo `pwd`
echo "==========================================="
ls -l
echo "==========================================="
echo
read -p "Ready to calculate?[y?]" judge
if [ "${judge}" == "y" ];then
  bsub<SUB;pwd;bjobs|awk 'END{ print }'
fi

exit 0
```



### 第二步

一旦判断出这是一个金属，那么，在这步优化之前要做一下SIGMA值的优化。

SIGMA值优化脚本

```shell
#!/bin/bash

echo "Welcome! Please type the range and precision of the SIGMA optimization:"
read -p "From:  " begin
read -p "to:  " end
read -p "Interval:  " prec
read -p "LMAXMIX:  " LMAXMIX_set
read -p "ENCUT: " ENCUT_set
read -p "ISMEAR:  " ISMEAR_set
echo
echo "======================================"
echo "Range:    ${begin} ~ ${prec} ~ ${end}"
echo "LMAXMIX:  ${LMAXMIX_set}"
echo "ENCUT:    ${ENCUT_set}"
echo "ISMEAR:   ${ISMEAR_set}"
echo "======================================"
echo
hom=`pwd`
rm -rf ${hom}/sigma_opt
mkdir ${hom}/sigma_opt
cp POSCAR POTCAR KPOINTS SUB ${hom}/sigma_opt
cd ${hom}/sigma_opt

hom=`pwd`
for i in $(seq ${begin} ${prec} ${end}) 
do
if [ -e "stop" ];then
  echo "The owner has kiiled the task!"
  echo "Have a good day!"
  exit 0
fi
direct_name=${hom}/sigma_${i}
rm -rf ${direct_name}
mkdir ${direct_name}
cd ${hom}
cp POSCAR POTCAR KPOINTS SUB ${direct_name}
cd ${direct_name}
cat > INCAR <<!
SYSTEM = SIGMA_test
ENCUT = ${ENCUT_set};LMAXMIX = ${LMAXMIX_set}
ISMEAR = ${ISMEAR_set} ; SIGMA = $i
NPAR = 1 ; PREC = Accurate
ISTART = 0 ; ICHARG = 2
LWAVE = .FALSE. ; LCHARG = .FALSE.
!

pwd
date
echo "SIGMA = $i eV"
bsub -K < SUB
tail -14 OUTCAR
judge=""
judge=`tail -14 OUTCAR|head -1|awk 'END{if($1=="General"){print 0}else{print 1}}'`
case ${judge} in
0)
  echo "The calculation in ${direct_name} successfully finished!"
  TS=`grep "EENTRO" OUTCAR | tail -1 | awk '{printf "%12.6f \n", $5 }'`
  echo $i $TS
  echo $i $TS >>${hom}/comment
;;
1)
  echo "Something goes wrong during the calculation in ${direct_name}."
;;
*)
  echo "Internal mistake!"
;;
esac
cd ${hom}
done
echo
echo "comment"
echo "==========================================="
cat ${hom}/comment 
echo "==========================================="
echo
echo "Successfully finished."
echo "Have a good day!"
exit 0
```



拷贝第一步的WAVECAR、CONTCAR到下一个工作目录

SUB

```shell
#BSUB -q e52680v3ib!
#BSUB -n 24
#BSUB -J 2_opt
#BSUB -o out
#BSUB -e err
module load ips/2019u5
mpirun vasp_std > log
```



KPOINTS

````shell
Automatic mesh
0
Gamma
2   2   2
0.0 0.0 0.0
````

INCAR

```shell
# Comments
SYSTEM = opt_2

# I/O
ISTART = 1
LWAVE = .FALSE.
LCHARG = .FALSE.

# Electronic Relaxation
ENCUT = 600
NELM = 100
PREC = Accurate
ISMEAR = 1 # -5 / 0 / 1
SIGMA = 0.2
EDIFF = 1E-7
LMAXMIX = 4
#LREAL = Auto
#ISYM = 0

# Ionic Relaxation
IBRION = 2
ISIF = 2
EDIFFG = 1E-4
NSW = 100

# Bands
NBANDS = 384
#LORBIT = 12

# Parallelization
NPAR = 1
```

记录一下计算的时间以及熵值，以及是否收敛（reach关键字）

P.S. 如果遇到如下错误，可以设置 **SYMPREC** **=** **1E-8** 或者 **ISYM=0** 解决.

```shell
Reciprocal lattice and k-lattice belong to different class of lattices. Often results are still useful...   60
```



## 静态计算

拷贝结构优化CONTCAR到下一个工作目录

SUB

```shell
#BSUB -q e52680v3ib!
#BSUB -n 24
#BSUB -J stat
#BSUB -o out
#BSUB -e err
module load ips/2019u5
mpirun vasp_std > log
```

KPOINTS

```shell 
Automatic mesh
0
Gamma
5   5   5
0.0 0.0 0.0
```

INCAR

```shell
# Comments
SYSTEM = stat

# I/O
ISTART = 0
LWAVE = .FALSE.
#LCHARG = .FALSE.

# Electronic Relaxation
ENCUT = 600
NELM = 100
PREC = Accurate
ISMEAR = -5
#SIGMA = 0.1
EDIFF = 1E-7
LMAXMIX = 4
#LREAL = Auto
#ISYM = 0

# Ionic Relaxation
#IBRION = 2
#ISIF = 2
#EDIFFG = 1E-4
#NSW = 100

# Bands
NBANDS = 192
LORBIT = 12

# Parallelization
NPAR = 1
```

静态计算脚本sta

```shell
#! /bin/bash

hom=`pwd`
dir="${hom}/stat"
mkdir ${dir}
cd ${hom}
cp -pai CONTCAR ${dir}/POSCAR
cp -pai POTCAR ${dir}
cd ${dir}
cat >INCAR<<!
# Comments
SYSTEM = stat

# I/O
ISTART = 0
LWAVE = .FALSE.
#LCHARG = .FALSE.

# Electronic Relaxation
ENCUT = 600
NELM = 100
PREC = Accurate
ISMEAR = -5
#SIGMA = 0.1
EDIFF = 1E-7
LMAXMIX = 4
#LREAL = Auto
#ISYM = 0

# Ionic Relaxation
#IBRION = 2
#ISIF = 2
#EDIFFG = 1E-4
#NSW = 100

# Bands
NBANDS = 192
LORBIT = 12

# Parallelization
NPAR = 1
!

cat >KPOINTS<<!
Automatic mesh
0
Gamma
5   5   5
0.0 0.0 0.0
!

cat >SUB<<!
#BSUB -q e52680v3ib!
#BSUB -n 24
#BSUB -J stat
#BSUB -o out
#BSUB -e err
module load ips/2019u5
mpirun vasp_std > log
!

echo
echo ${dir}
echo "==========================================="
ls -l
echo "==========================================="
echo

bsub<SUB
pwd;bjobs|awk 'END{ print }'

exit 0
```

## 能带计算

221超胞的INCAR文件

```shell
# Comments
SYSTEM = dos

# I/O
ISTART = 0
ICHARG = 11
LWAVE = .FALSE.
LCHARG = .FALSE.

# Electronic Relaxation
ENCUT = 600
NELM = 100
PREC = Accurate
ISMEAR = 0
#SIGMA = 0.1
EDIFF = 1E-7
LMAXMIX = 4
#LREAL = Auto
#ISYM = 0

# Ionic Relaxation
#IBRION = 2
#ISIF = 2
#EDIFFG = 1E-4
#NSW = 100

# Bands
NBANDS = 192
LORBIT = 12

# Parallelization
NPAR = 1
```

KPOINTS文件???由Vaspkit生成

```shell
k-points along high symmetry lines 
120           ! intersections 24*5
Line-mode 
reciprocal 
0.0 0.0 0.0 !G
0.5 0.5 0.0 !X

0.5 0.5 0.0 !X
0.5 0.5 0.25 !P

0.5 0.5 0.25 !P
0.5 0.0 0.0 !N

0.5 0.0 0.0 !N
0.0 0.0 0.0 !G

0.0 0.0 0.0 !G
0.0 0.0 0.5 !Z
```

SUB

```shell
#BSUB -q e52680v3ib!
#BSUB -n 24
#BSUB -J dos
#BSUB -o out
#BSUB -e err
module load ips/2019u5
mpirun vasp_std > log
```

## 态密度计算

221超胞的INCAR文件

```shell
# Comments
SYSTEM = dos

# I/O
ISTART = 0
ICHARG = 11
LWAVE = .FALSE.
LCHARG = .FALSE.

# Electronic Relaxation
ENCUT = 600
NELM = 100
PREC = Accurate
ISMEAR = -5
EDIFF = 1E-7
LMAXMIX = 4

# DOS: self_EF =
NEDOS = 1500  # 301-Default  2000~3000
EMIN =    # EF-1
EMAX =    # EF+1

# Bands
NBANDS = 192
LORBIT = 12  # 11-DOSCAR+lm-decomposed PROCAR 12-DOSCAR+lm-decomposed PROCAR+phase factors

# Parallelization
NPAR = 1
```

KPOINTS

```shell
Automatic mesh
0
Gamma
10  10   10
0.0 0.0 0.0
```

SUB

```shell
#BSUB -q e52680v3ib!
#BSUB -n 24
#BSUB -J dos
#BSUB -o out
#BSUB -e err
module load ips/2019u5
mpirun vasp_std > log
```



判断两个POSCAR文件是否等价

```shell
#! /bin/bash

hom=`pwd`
declare -a file
for i in 1 2
do
  read -p "file $i: " file[$i]
  if [ ! -e "${hom}/${file[$i]}" ];then
    echo "There is no ${file[$i]}."
    echo "Exit."
    echo "Have a good day!"
    exit 1
fi
done
dir="${hom}/compare"
a=1
while [ -e ${dir} ]
do
  dir="${hom}/compare_${a}"
  let a++
done
mkdir ${dir}
echo
cp -pai ${file[1]} ${dir}/POSCAR_1
cp -pai ${file[2]} ${dir}/POSCAR_2
file1_p=`printf "%16s" ${file[1]}`
file2_p=`printf "%16s" ${file[2]}`
echo
echo "compare"
echo "==========================================="
echo "${file1_p} =======>  POSCAR_1"
echo "${file1_p} =======>  POSCAR_2"
echo "==========================================="
echo
cd ${dir}

for i in 1 2
do
  phonopy -c POSCAR_${i} --symmetry>sym_${i}
  mv PPOSCAR POSCAR_g_${i}
  phonopy -c POSCAR_g_${i} --symmetry>sym_${i}
done

file1="sym_1"
file2="sym_2"

lines1=`cat ${file1} | wc -l`
lines2=`cat ${file2} | wc -l`
if (( ${lines1} != ${lines2} ));then
  echo
  echo "POSCAR files are not equivalent."
  echo
  echo "Finished."
  echo "Have a good day!"
  exit 0
fi

for ((i=1;i<=${lines1};i++))
do
  line1=`awk 'NR=="'$i'"{print $0}' ${file1}`
  line2=`awk 'NR=="'$i'"{print $0}' ${file2}`
  if [[ "${line1}" != "${line2}" ]];then
    echo
    echo "POSCAR files may not be equivalent."
    echo
    echo "Finished."
    echo "Have a good day!"
    exit 0
  fi
done
echo
echo "POSCAR files are equivalent."
echo
echo "Finished."
echo "Have a good day!"

exit 0
```

判断多个POSCAR文件之间是否等价

```shell
#! /bin/bash

hom=`pwd`
read -p "Which directory: " temp
hom="${hom}/${temp}"
declare -a file
i=1
if [ -e "${hom}/order" ];then
  rm ${hom}/order
fi
if [ -e "${hom}/equi.out" ];then
  rm ${hom}/equi.out
fi
echo
for temp in `ls ${hom}`
do
  file[$i]="${temp}"
  printf "%4d%16s\n" $i ${temp}>>${hom}/order
  let i++
done
num=${#file[*]}
echo $num
for ((i=1;i<${num};i++))
do
  for ((j=${i}+1;j<=${num};j++))
  do
    let order=(i-1)*num+j
    judge[${order}]=0
    dir="${hom}/compare_${i}${j}"
    mkdir ${dir}
    echo
    cp -pai ${hom}/${file[i]} ${dir}/POSCAR_${i}
    cp -pai ${hom}/${file[j]} ${dir}/POSCAR_${j}
    file1_p=`printf "%16s" ${file[i]}`
    file2_p=`printf "%16s" ${file[j]}`
    echo
    echo "compare_${i}${j}"
    echo "==========================================="
    echo "${file1_p} =======>  POSCAR_${i}"
    echo "${file2_p} =======>  POSCAR_${j}"
    echo "==========================================="
    echo

    cd ${dir}

    for k in $i $j
    do
      phonopy -c POSCAR_${k} --symmetry>sym_${k}
      mv PPOSCAR POSCAR_g_${k}
      phonopy -c POSCAR_g_${k} --symmetry>sym_${k}
    done

    file1="POSCAR_g_${i}"
    file2="POSCAR_g_${j}"

    file1="sym_${i}"
    file2="sym_${j}"

    lines1=`cat ${file1} | wc -l`
    lines2=`cat ${file2} | wc -l`
    if (( ${lines1} != ${lines2} ));then
      judge[${order}]=1
      continue
    fi

    for ((k=1;k<=${lines1};k++))
    do
      line1=`awk 'NR=="'$k'"{print $0}' ${file1}`
      line2=`awk 'NR=="'$k'"{print $0}' ${file2}`
      if [[ "${line1}" != "${line2}" ]];then
        judge[${order}]=1
        continue
      fi
    done
  cd ${hom}
  done
done

echo
echo "Equivalent file"
echo "==========================================="
for ((i=1;i<${num};i++))
do
  for ((j=${i}+1;j<=${num};j++))
  do
    let order=(i-1)*num+j
    if [ "${judge[${order}]}" == "0" ];then
      printf "%16s%16s\n" ${file[$i]} ${file[$j]}| tee -a ${hom}/equi.out
    fi
  done
done
echo "==========================================="
echo
echo "Finished."
echo "Have a good day!"

exit 0
```



改进版本（烂尾楼）

```shell
#! /bin/bash

hom=`pwd`
read -p "Which directory: " temp
hom="${hom}/${temp}"
declare -a file
i=1
if [ -e "${hom}/order" ];then
  rm ${hom}/order
fi
if [ -e "${hom}/equi.out" ];then
  rm ${hom}/equi.out
fi
echo
for temp in `ls ${hom}`
do
  file[$i]="${temp}"
  printf "%4d%16s\n" $i ${temp}>>${hom}/order
  let i++
done
num=${#file[*]}
echo $num
for ((i=1;i<${num};i++))
do
  for ((j=${i}+1;j<=${num};j++))
  do
    let order=(i-1)*num+j
    judge[${order}]=0
    dir="${hom}/compare_${i}${j}"
    mkdir ${dir}
    echo
    cp -pai ${hom}/${file[i]} ${dir}/POSCAR_${i}
    cp -pai ${hom}/${file[j]} ${dir}/POSCAR_${j}
    file1_p=`printf "%16s" ${file[i]}`
    file2_p=`printf "%16s" ${file[j]}`
    echo
    echo "compare_${i}${j}"
    echo "==========================================="
    echo "${file1_p} =======>  POSCAR_${i}"
    echo "${file2_p} =======>  POSCAR_${j}"
    echo "==========================================="
    echo

    cd ${dir}

    for k in $i $j
    do
      phonopy -c POSCAR_${k} --symmetry>sym_${k}
    done
    
    # 先写i的，把所有信息写进一个文件中
    k=1
    for temp in `awk 'NR==6{print $0}' POSCAR_${i}`
    do
      atom_type[k]=${temp}  # 原子种类
      let k++
    done
    k=1
    for temp in `awk 'NR==7{print $0}' POSCAR_${i}`
    do
      atom_num[k]=${temp}  # 原子数目
      let total_num=total_num+temp
      let k++
    done
    t=1
    group_num=`awk 'NR==3{print $2}' sym_${t}` #空间群群号
    mark=`awk '$1=="atom_mapping:"{print NR}' sym_${t}`
    let upp=total_num+mark
    k=1
    for ((m=${mark};m<=${upp};m++))
    do
      map[${k}]=`awk -v m=${m} 'NR==m{print $2}' sym_${t}`
      let k++
    done
    r=1
    for ((k=1;k<=${total_num};k++))
    do
      q=0
      for ((m=1;m<=${total_num};m++))
      do
        if (( ${map[${m}]} == $k ));then
          let q++
        fi
      done
      if (( $q != 0 ));then
        atom_irre[$r]=$k
        equiv_num[$r]=${q}
        Wyckoff_mark[$r]=`grep "Wyckoff" sym_${t}|awk 'NR=="'${r}'"{print $2}'`
        let r++
      fi
    
    eval atom${i}_irre[$r]=${atom_irre[$r]}
    
    done
    
  cd ${hom}
  done
done

echo
echo "Equivalent file"
echo "==========================================="
for ((i=1;i<${num};i++))
do
  for ((j=${i}+1;j<=${num};j++))
  do
    let order=(i-1)*num+j
    if [ "${judge[${order}]}" == "0" ];then
      printf "%16s%16s\n" ${file[$i]} ${file[$j]}| tee -a ${hom}/equi.out
    fi
  done
done
echo "==========================================="
echo
echo "Finished."
echo "Have a good day!"

exit 0
```

最终

```shell
#! /bin/bash
# This is a script to judge the equivalence among POSCAR files provided in the selected folder.

##################################################################
########################  function part   ########################
##################################################################
function same(){
  fi=1
  for tempo in $@
  do
    case $fi in
    1)
      num1=${tempo}
      let fi++
    ;;
    2)
      num2=${tempo}
      let fi++
    ;;
    3)
      prec=1
      for ((fj=1;fj<=${tempo};fj++))
      do
        prec="${prec}0"
      done   
      let fi++
    ;;
    *)
      break
    ;;
    esac
  done
  abs=`awk -v p=${prec} -v n1=${num1} -v n2=${num2} 'BEGIN{a=p*(n1-n2);print a}'` 
  abs=`awk -v abs=${abs} 'BEGIN{print (abs<0)?(-abs):(abs)}'`
  judge=`awk -v abs=${abs} 'BEGIN{print (abs<1)?0:1}'`
  return ${judge}
}
##################################################################
##################################################################

##################################################################
##########################  Input part   #########################
##################################################################
echo
echo "Welcome!"
echo "This is a script to judge the equivalence among POSCAR files provided in the selected folder."
echo
read -p "Which directory: " temp
hom=`pwd`
work="${hom}/workout"
rm -rf ${work}
mkdir ${work}
hom="${hom}/${temp}"
echo
rm -rf ${hom}/compare_*
pre=7
read -p "Precision: " temp
if [ "${temp}" != "" ];then
  pre=${temp}
fi
touch ${work}/log
echo |tee -a ${work}/log
echo "Precision: ${pre}"|tee -a ${work}/log
echo|tee -a ${work}/log
echo "File orders"|tee -a ${work}/log
echo "==========================================="|tee -a ${work}/log
echo|tee -a ${work}/log
declare -a file
i=1
for temp in `ls ${hom}`
do
  file[$i]="${temp}"
  printf "%4d%16s\n" $i ${temp}>>${work}/order.out
  let i++
done
cat ${work}/order.out >> ${work}/log
echo "==========================================="|tee -a ${work}/log
num=${#file[*]}
##################################################################
##################################################################

##################################################################
######################   Data process part   #####################
##################################################################
for ((i=1;i<${num};i++))
do
  let ip=i+1
  for ((j=${ip};j<=${num};j++))
  do
    let order=(i-1)*num+j
    judge[${order}]=0
    dir="${hom}/compare_${i}${j}"
    mkdir ${dir}
    cp -pai ${hom}/${file[$i]} ${dir}/POSCAR_${i}
    cp -pai ${hom}/${file[$j]} ${dir}/POSCAR_${j}
    file1_p=`printf "%16s" ${file[$i]}`
    file2_p=`printf "%16s" ${file[$j]}`
    echo|tee -a ${work}/log
    echo "compare_${i}${j}"|tee -a ${work}/log
    echo "==========================================="||tee -a ${work}/log
    echo "${file1_p} =======>  POSCAR_${i}"|tee -a ${work}/log
    echo "${file2_p} =======>  POSCAR_${j}"|tee -a ${work}/log
    echo "==========================================="|tee -a ${work}/log
    echo|tee -a ${work}/log
    cd ${dir}
    for k in $i $j
    do
      phonopy -c POSCAR_${k} --symmetry>sym_${k}
      mv BPOSCAR POSCAR_g_${k}
    done
    file1="POSCAR_g_${i}"
    file2="POSCAR_g_${j}"
    for k in 3 4 5
    do
      for l in 1 2 3
      do
        num1=`awk -v l=${l} 'NR=="'$k'"{print $l}' ${file1}`
        num2=`awk -v l=${l} 'NR=="'$k'"{print $l}' ${file2}`
        same ${num1} ${num2} ${pre}
        temp=$?
        let judge[${order}]=judge[${order}]+temp
      done
    done
    temp=`awk 'ARGIND==1{a=0;if(NR==7){a=$1+$2+$3+$4}}ARGIND==2{b=0;if(NR==7){b=$1+$2+$3+$4}}END{print (a!=b)?1:0}' ${file1} ${file2}`
    let judge[${order}]=judge[${order}]+temp
    if [ "${judge[${order}]}" != "0" ];then
      continue
    fi
    num_true=`awk 'END{print NR}' ${file1}`
    for ((k=9;k<=${num_true};k++))
    do
      for l in 1 2 3
      do
        num1=`awk -v l=${l} 'NR=="'$k'"{print $l}' ${file1}`
        num2=`awk -v l=${l} 'NR=="'$k'"{print $l}' ${file2}`
        same ${num1} ${num2} ${pre}
        temp=$?
        let judge[${order}]=judge[${order}]+temp
        if [ "${judge[${order}]}" != "0" ];then
          break
        fi
      done
    done
  cd ${hom}
  done
done
##################################################################
##################################################################

##################################################################
##########################   Output part   #######################
##################################################################
echo|tee -a ${work}/log
echo "Equivalent file"|tee -a ${work}/log
echo "==========================================="|tee -a ${work}/log
for ((i=1;i<${num};i++))
do
  for ((j=${i}+1;j<=${num};j++))
  do
    let order=(i-1)*num+j
    if [ "${judge[${order}]}" == "0" ];then
      printf "%16s%16s\n" ${file[$i]} ${file[$j]}| tee -a ${work}/equi.out
    fi
  done
done
cat ${work}/equi.out >>${work}/log
echo "==========================================="|tee -a ${work}/log
echo|tee -a ${work}/log
echo "All outputs have been saved to the directory:"|tee -a ${work}/log
echo ${work}|tee -a ${work}/log
echo|tee -a ${work}/log
echo "Finished."|tee -a ${work}/log
echo "Have a good day!"|tee -a ${work}/log
echo
##################################################################
##################################################################
exit 0
```



sym文件读取

```shell
 #! /bin/bash
   
   read i
   # 先写i的，把所有信息写进一个文件中
    k=0
    for temp in `awk 'NR==6{print $0}' POSCAR_${i}`
    do
      let k++
      atom_type[$k]=${temp}  # 原子种类
    done
    k=0
    total_num=0
    for temp in `awk 'NR==7{print $0}' POSCAR_${i}`
    do
      let k++
      atom_num[$k]=${temp}  # 原子数目
      let total_num=total_num+temp
    done
    group_num=`awk 'NR==3{print $2}' sym_${i}` # 空间群群号
    mark=`awk '$1=="atom_mapping:"{print NR}' sym_${i}` # 定位sym文件中的等价位置信息
    let mark++
    let upp=total_num+mark-1
    k=0
    for ((m=mark;m<=upp;m++))
    do
      let k++
      map[${k}]=`awk 'NR=="'$m'"{print $2}' sym_${i}` 
    done
    r=0
    for ((k=1;k<=total_num;k++))
    do
      q=0
      for ((m=1;m<=total_num;m++))
      do
        if (( map[m] == k ));then
          let q++
        fi
      done
      if (( q != 0 ));then
        let r++
        irre[$r]=$k  # 每个不等价原子在POSCAR中的序号
        equiv_num[$r]=${q}  # 每个不等价原子重数
        Wyckoff_mark[$r]=`grep "Wyckoff" sym_${i}|awk 'NR=="'${r}'"{print $2}'|sed $'s/\'//g'` # 每个不等价原子的Wyckoff标记
      fi
    done
    limit_bot=1
    for ((type=1;type<=${#atom_num[*]};type++))
    do
      let limit_up=limit_bot+atom_num[${type}]-1
      temp=0;order=0;
      for ((k=1;k<=r;k++))
      do
        if (( irre[$k] >= limit_bot && irre[$k] <= limit_up ));then
          let order++;let temp++;
          eval atom${type}_irre['$'{order}]='$'{irre['$'k]}
          eval atom${type}_equiv_num['$'{order}]='$'{equiv_num['$'k]}
          eval atom${type}_Wyckoff_mark['$'{order}]='$'{Wyckoff_mark['$'k]}
        fi
      done
      type_num[${type}]=${temp}  # 每类原子的个数
      let limit_bot=limit_up+1
    done
    for ((type=1;type<=${#atom_num[*]};type++))
    do
      eval declare -a atom${type}_order
      for ((k=1;k<=type_num[${type}];k++))
      do 
        eval printf "%3s%3d%3d%4s'\n'" '$'{atom_type['$'type]} '$'{atom${type}_irre['$'k]} '$'{atom${type}_equiv_num['$'k]} '$'{atom${type}_Wyckoff_mark['$'k]} 
        eval atom${type}_order['$'k]='$'k
      done
    done

# 排序，先按照Wyckoff记号，a b c d ..., 然后按照原子的数目从低到高来排列，如果遇到相等的情况，按照原先位置排列，也就是大于才更换位置
    for ((type=1;type<=${#atom_num[*]};type++))
    do
      eval declare -a atom${type}_order
      for ((k1=1;k1<=type_num[${type}];k1++))
      do 
        eval w1='`'printf "%d" '\'"'""$"{atom${type}_Wyckoff_mark[$k1]}'`' 
        for ((k2=1;k2<=type_num[${type}];k2++))
        do
          eval w2='`'printf "%d" '\'"'""$"{atom${type}_Wyckoff_mark[$k2]}'`' 
          if (( w2 < w1 ));then
            eval o1='$'{atom${type}_order['$'k1]}
            eval o2='$'{atom${type}_order['$'k2]}
            if (( o2 > o1 ));then
              eval atom${type}_order['$'k1]='$'o2
              eval atom${type}_order['$'k2]='$'o1 
            fi
          elif (( w2 == w1 ));then
            eval x1='`'printf "%d" "$"{atom${type}_equiv_num[$k1]}'`' 
            eval x2='`'printf "%d" "$"{atom${type}_equiv_num[$k2]}'`'
            if (( x2 < x1 ));then
              eval o1='$'{atom${type}_order['$'k1]}
              eval o2='$'{atom${type}_order['$'k2]}
              if (( o2 > o1 ));then
                eval atom${type}_order['$'k1]='$'o2
                eval atom${type}_order['$'k2]='$'o1 
              fi
            fi
          fi
        done
      done
    done

    for ((type=1;type<=${#atom_num[*]};type++))
    do
      for ((k=1;k<=type_num[${type}];k++))
      do 
        eval printf "%3s%3s%3d%3d%4s'\n'" '$'{atom${type}_order['$'k]} '$'{atom_type['$'type]} '$'{atom${type}_irre['$'k]} '$'{atom${type}_equiv_num['$'k]} '$'{atom${type}_Wyckoff_mark['$'k]} 
      done
    done

exit 0
```



最终版1.0

```shell
#! /bin/bash
# This is a script to judge the equivalence among POSCAR files provided in the selected folder.

##################################################################
########################  function part   ########################
##################################################################
function same(){
  fi=1
  for tempo in $@
  do
    case $fi in
    1)
      n1=${tempo}
      let fi++
    ;;
    2)
      n2=${tempo}
      let fi++
    ;;
    3)
      prec=1
      for ((fj=1;fj<=${tempo};fj++))
      do
        prec="${prec}0"
      done   
      let fi++
    ;;
    *)
      break
    ;;
    esac
  done
  abs=`awk -v p=${prec} -v n1=${n1} -v n2=${n2} 'BEGIN{a=p*(n1-n2);print a}'` 
  abs=`awk -v abs=${abs} 'BEGIN{print (abs<0)?(-abs):(abs)}'`
  jud=`awk -v abs=${abs} 'BEGIN{print (abs<1)?0:1}'`
  return ${jud}
}
##################################################################
##################################################################

##################################################################
##########################  Input part   #########################
##################################################################
echo
echo "Welcome!"
echo "This is a script to judge the equivalence among POSCAR files provided in the selected folder."
echo
read -p "Which directory: " temp
hom=`pwd`
work="${hom}/workout"
rm -rf ${work}
mkdir ${work}
hom="${hom}/${temp}"
echo
rm -rf ${hom}/compare_*
pre=7
read -p "Precision: " temp
if [ "${temp}" != "" ];then
  pre=${temp}
fi
touch ${work}/log
touch ${work}/order.out
touch ${work}/equi.out
echo |tee -a ${work}/log
echo "Precision: ${pre}"|tee -a ${work}/log
echo|tee -a ${work}/log
echo "File orders"|tee -a ${work}/log
echo "==========================================="|tee -a ${work}/log
echo|tee -a ${work}/log
declare -a file
i=1
for temp in `ls ${hom}`
do
  file[$i]="${temp}"
  printf "%4d%16s\n" $i ${temp}| tee -a ${work}/order.out
  let i++
done
cat ${work}/order.out >> ${work}/log
echo "==========================================="|tee -a ${work}/log
num=${#file[*]}
##################################################################
##################################################################

##################################################################
######################   Data process part   #####################
##################################################################
for ((i=1;i<num;i++))
do
  for ((j=i+1;j<=num;j++))
  do
    let order=(i-1)*num+j
    judge[${order}]=0
    dir="${hom}/compare_${i}${j}"
    mkdir ${dir}
    cp -pai ${hom}/${file[$i]} ${dir}/POSCAR_${i}
    cp -pai ${hom}/${file[$j]} ${dir}/POSCAR_${j}
    file1_p=`printf "%16s" ${file[$i]}`
    file2_p=`printf "%16s" ${file[$j]}`
    echo|tee -a ${work}/log
    echo "compare_${i}${j}"|tee -a ${work}/log
    echo "==========================================="|tee -a ${work}/log
    echo "${file1_p} =======>  POSCAR_${i}"|tee -a ${work}/log
    echo "${file2_p} =======>  POSCAR_${j}"|tee -a ${work}/log
    echo "==========================================="|tee -a ${work}/log
    echo|tee -a ${work}/log
    cd ${dir}
    for k in $i $j
    do
      phonopy -c POSCAR_${k} --symmetry>sym_${k}
      mv PPOSCAR POSCAR_g_${k}
    done
    file1="POSCAR_g_${i}"
    file2="POSCAR_g_${j}"
    # 首先比较一下phonopy转化生成的原胞的基矢是不是一致。phonopy生成的原胞估计是考虑了点群，但没有考虑滑移
    for k in 3 4 5
    do
      for l in 1 2 3
      do
        num1=`awk -v l=${l} 'NR=="'$k'"{print $l}' ${file1}`
        num2=`awk -v l=${l} 'NR=="'$k'"{print $l}' ${file2}`
        same ${num1} ${num2} ${pre}
        temp=$?
        let judge[${order}]=judge[${order}]+temp
      done
    done
    # 比较一下原子总数是否一致
    temp=`awk 'ARGIND==1{a=0;if(NR==7){a=$1+$2+$3+$4}}ARGIND==2{b=0;if(NR==7){b=$1+$2+$3+$4}}END{print (a!=b)?1:0}' ${file1} ${file2}`
    let judge[${order}]=judge[${order}]+temp
    if (( judge[${order}] != 0 ));then
      continue
    fi

    for ii in $i $j
    do
       # 先写i后写j的，把两个晶胞的对称性信息分别写进两个文件中
      k=0
      for temp in `awk 'NR==6{print $0}' POSCAR_g_${ii}`
      do
        let k++
        atom_type[$k]=${temp}  # 原子种类
      done
      k=0
      total_num=0
      for temp in `awk 'NR==7{print $0}' POSCAR_g_${ii}`
      do
        let k++
        atom_num[$k]=${temp}  # 原子数目
        let total_num=total_num+temp
      done
      group_num=`awk 'NR==3{print $2}' sym_${ii}` # 空间群群号
      mark=`awk '$1=="atom_mapping:"{print NR}' sym_${ii}` # 定位sym文件中的等价位置信息
      let mark++
      let upp=total_num+mark-1
      k=0
      for ((m=mark;m<=upp;m++))
      do
        let k++
        map[${k}]=`awk 'NR=="'$m'"{print $2}' sym_${ii}` 
      done
      r=0
      for ((k=1;k<=total_num;k++))
      do
        q=0
        for ((m=1;m<=total_num;m++))
        do
          if (( map[m] == k ));then
            let q++
          fi
        done
        if (( q != 0 ));then
          let r++
          irre[$r]=$k  # 每个不等价原子在POSCAR中的序号
          equiv_num[$r]=${q}  # 每个不等价原子重数
          Wyckoff_mark[$r]=`grep "Wyckoff" sym_${ii}|awk 'NR=="'${r}'"{print $2}'|sed $'s/\'//g'` # 每个不等价原子的Wyckoff标记
        fi
      done
      limit_bot=1
      # 把每种原子信息拆开放在不同数组
      for ((type=1;type<=${#atom_num[*]};type++))
      do
        let limit_up=limit_bot+atom_num[${type}]-1
        temp=0;orde=0;
        for ((k=1;k<=r;k++))
        do
          if (( irre[$k] >= limit_bot && irre[$k] <= limit_up ));then
          let orde++;let temp++;
          eval atom${type}_irre['$'{orde}]='$'{irre['$'k]}
          eval atom${type}_equiv_num['$'{orde}]='$'{equiv_num['$'k]}
          eval atom${type}_Wyckoff_mark['$'{orde}]='$'{Wyckoff_mark['$'k]}
          fi
        done
        type_num[${type}]=${temp}  # 每类原子的个数
        let limit_bot=limit_up+1
      done
      for ((type=1;type<=${#atom_num[*]};type++))
      do
        eval declare -a atom${type}_order
        for ((k=1;k<=type_num[${type}];k++))
        do 
          eval atom${type}_order['$'k]='$'k
        done
      done

# 排序：先按照Wyckoff记号，a b c d ..., 然后按照每个Wyckoff位置含有的原子的数目从低到高排列
      for ((type=1;type<=${#atom_num[*]};type++))
      do
        for ((k1=1;k1<=type_num[${type}]-1;k1++))
        do 
          eval w1='`'printf "%d" '\'"'""$"{atom${type}_Wyckoff_mark[$k1]}'`' 
          for ((k2=k1;k2<=type_num[${type}];k2++))
          do
            eval w2='`'printf "%d" '\'"'""$"{atom${type}_Wyckoff_mark[$k2]}'`' 
            if (( w2 < w1 ));then
              eval temp='$'{atom${type}_irre['$'k1]}
              eval atom${type}_irre['$'k1]='$'{atom${type}_irre['$'k2]}
              eval atom${type}_irre['$'k2]='$'{temp}            
              eval temp='$'{atom${type}_equiv_num['$'k1]}
              eval atom${type}_equiv_num['$'k1]='$'{atom${type}_equiv_num['$'k2]}
              eval atom${type}_equiv_num['$'k2]='$'{temp}       
              eval temp='$'{atom${type}_Wyckoff_mark['$'k1]}
              eval atom${type}_Wyckoff_mark['$'k1]='$'{atom${type}_Wyckoff_mark['$'k2]}
              eval atom${type}_Wyckoff_mark['$'k2]='$'{temp}
            elif (( w2 == w1 ));then
              eval x1='$'{atom${type}_equiv_num[$k1]}
              eval x2='$'{atom${type}_equiv_num[$k2]}
              if (( x2 < x1 ));then
                eval temp='$'{atom${type}_irre['$'k1]}
                eval atom${type}_irre['$'k1]='$'{atom${type}_irre['$'k2]}
                eval atom${type}_irre['$'k2]='$'{temp}            
                eval temp='$'{atom${type}_equiv_num['$'k1]}
                eval atom${type}_equiv_num['$'k1]='$'{atom${type}_equiv_num['$'k2]}
                eval atom${type}_equiv_num['$'k2]='$'{temp}       
                eval temp='$'{atom${type}_Wyckoff_mark['$'k1]}
                eval atom${type}_Wyckoff_mark['$'k1]='$'{atom${type}_Wyckoff_mark['$'k2]}
                eval atom${type}_Wyckoff_mark['$'k2]='$'{temp}
#              elif (( x2 == x1 ));then
#                eval y1='$'{atom${type}_irre[$k1]}
#                eval y2='$'{atom${type}_irre[$k2]}
#                if (( y2 < y1 ));then
#                  eval temp='$'{atom${type}_irre['$'k1]}
#                  eval atom${type}_irre['$'k1]='$'{atom${type}_irre['$'k2]}
#                  eval atom${type}_irre['$'k2]='$'{temp}            
#                  eval temp='$'{atom${type}_equiv_num['$'k1]}
#                  eval atom${type}_equiv_num['$'k1]='$'{atom${type}_equiv_num['$'k2]}
#                  eval atom${type}_equiv_num['$'k2]='$'{temp}       
#                  eval temp='$'{atom${type}_Wyckoff_mark['$'k1]}
#                  eval atom${type}_Wyckoff_mark['$'k1]='$'{atom${type}_Wyckoff_mark['$'k2]}
#                  eval atom${type}_Wyckoff_mark['$'k2]='$'{temp}
#                fi
              fi
            fi
          done
        done
      done  
      touch ${dir}/sym_log_${ii}
      echo ${group_num}>>${dir}/sym_log_${ii}
      for ((type=1;type<=${#atom_num[*]};type++))
      do
        for ((k=1;k<=type_num[${type}];k++))
        do      
          eval printf "%3s%4s%3s%3d'\n'" '$'{atom${type}_order['$'k]} '$'{atom_type['$'{type}]} '$'{atom${type}_Wyckoff_mark['$'k]} '$'{atom${type}_equiv_num['$'k]} |tee -a ${dir}/sym_log_${ii}
        done
      done
    done
     # 比较生成的两个sym_log文件是不是一样的
    file1=${dir}/sym_log_${i}
    file2=${dir}/sym_log_${j}
    # 先比较行数
    lines1=`cat ${file1} | wc -l`
    lines2=`cat ${file2} | wc -l`
    if (( lines1 != lines2 ));then
      let judge[${order}]++
      continue
    fi
    # 再比较细节
    for ((k=1;k<=lines1;k++))
    do  
      line1=`awk 'NR=="'$k'"{print $0}' ${file1}`
      line2=`awk 'NR=="'$k'"{print $0}' ${file2}` 
      if [[ "${line1}" != "${line2}" ]];then
        let judge[${order}]++
        break 
      fi    
    done
  cd ${hom}
  done
done
##################################################################
##################################################################

##################################################################
##########################   Output part   #######################
##################################################################
echo|tee -a ${work}/log
echo "Equivalent file"|tee -a ${work}/log
echo "==========================================="|tee -a ${work}/log
for ((i=1;i<num;i++))
do
  for ((j=i+1;j<=num;j++))
  do
    let order=(i-1)*num+j
    if [ "${judge[${order}]}" == "0" ];then
      printf "%16s%16s\n" ${file[$i]} ${file[$j]}| tee -a ${work}/equi.out
    fi
  done
done
cat ${work}/equi.out >>${work}/log
echo "==========================================="|tee -a ${work}/log
echo|tee -a ${work}/log
echo "All outputs have been saved to the directory:"|tee -a ${work}/log
echo ${work}|tee -a ${work}/log
echo|tee -a ${work}/log
echo "Finished."|tee -a ${work}/log
echo "Have a good day!"|tee -a ${work}/log
echo
##################################################################
##################################################################
exit 0
```

