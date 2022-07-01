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
# poseq
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
hom=`pwd`
echo
echo "Welcome!"
echo "This is a script to judge the equivalence among POSCAR files provided in the selected folder."
echo
if [[ "$1" != "" ]];then
  temp=$1
else
  echo
  echo "Current work directory :" ${hom}
  echo "==========================================="
  ls
  echo "==========================================="
  echo
  read -p "Which directory: " temp  #输入当前工作目录下的POSCAR文件所在目录
fi
if [[ ! -e "${hom}/${temp}" ]];then
  echo "No such file!"
  echo
  echo "Quit"
  echo "Have a good day!"
  exit 1
fi
pre=7
if [[ "$2" != "" ]];then
  pre=$2
else
  read -p "Precision: " temp
  if [ "${temp}" != "" ];then
    pre=${temp}
  fi
fi
work="${hom}/workout"
rm -rf ${work}
mkdir ${work}
hom="${hom}/${temp}"
echo
rm -rf ${hom}/compare_*
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
      for temp in `awk 'NR==6{print $0}' POSCAR_${ii}`
      do
        let k++
        atom_type[$k]=${temp}  # 原子种类
      done
      k=0
      total_num=0
      for temp in `awk 'NR==7{print $0}' POSCAR_${ii}`
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
          for ((k2=k1+1;k2<=type_num[${type}];k2++))
          do
            eval w1='`'printf "%d" '\'"'""$"{atom${type}_Wyckoff_mark[$k1]}'`'
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
              elif (( x2 == x1 ));then
                eval y1='$'{atom${type}_irre[$k1]}
                eval y2='$'{atom${type}_irre[$k2]}
                if (( y2 < y1 ));then
                  eval temp='$'{atom${type}_irre['$'k1]}
                  eval atom${type}_irre['$'k1]='$'{atom${type}_irre['$'k2]}
                  eval atom${type}_irre['$'k2]='$'{temp}            
                  eval temp='$'{atom${type}_equiv_num['$'k1]}
                  eval atom${type}_equiv_num['$'k1]='$'{atom${type}_equiv_num['$'k2]}
                  eval atom${type}_equiv_num['$'k2]='$'{temp}       
                  eval temp='$'{atom${type}_Wyckoff_mark['$'k1]}
                  eval atom${type}_Wyckoff_mark['$'k1]='$'{atom${type}_Wyckoff_mark['$'k2]}
                  eval atom${type}_Wyckoff_mark['$'k2]='$'{temp}
                fi
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



判断两个POSCAR文件是不是一模一样

```shell
#! /bin/bash

##################################################################
########################  function part   ########################
##################################################################
function same(){
  abs=`awk -v p=$3 -v n1=$1 -v n2=$2 'BEGIN{a=p*(n1-n2);print a}'` 
  abs=`awk -v abs=${abs} 'BEGIN{print (abs<0)?(-abs):(abs)}'`
  judge=`awk -v abs=${abs} 'BEGIN{print (abs<1)?0:1}'`
  return ${judge}
}
##################################################################
##################################################################
##################################################################

pre=7
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
  echo $i
  let ip=i+1
  for ((j=${ip};j<=${num};j++))
  do
    echo $j
    let order=(i-1)*num+j
    judge[${order}]=0
    dir="${hom}/compare_${i}${j}"
    mkdir ${dir}
    echo
    cp -pai ${hom}/${file[$i]} ${dir}/POSCAR_${i}
    cp -pai ${hom}/${file[$j]} ${dir}/POSCAR_${j}
    file1_p=`printf "%16s" ${file[$i]}`
    file2_p=`printf "%16s" ${file[$j]}`
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

```shell
#! /bin/bash
# possame
# This is a script to judge the smae POSCAR files provided in the selected folder.

##################################################################
########################  function part   ########################
##################################################################
function same(){
  prec=1
  for ((fj=1;fj<=$3;fj++))
  do
    prec="${prec}0"
  done   
  abs=`awk -v p=${prec} -v n1=$1 -v n2=$2 'BEGIN{a=p*(n1-n2);print a}'` 
  abs=`awk -v abs=${abs} 'BEGIN{print (abs<0)?(-abs):(abs)}'`
  judge_1=`awk -v abs=${abs} 'BEGIN{print (abs<1)?0:1}'`
  abs=`awk -v p=${prec} -v abs=${abs} 'BEGIN{a=p-abs;print a}'`
  abs=`awk -v abs=${abs} 'BEGIN{print (abs<0)?(-abs):(abs)}'`  
  judge_2=`awk -v abs=${abs} 'BEGIN{print (abs<1)?0:1}'`
  let judge=(judge_1+judge_2)/2
  return ${judge}
}
##################################################################
##################################################################

##################################################################
##########################  Input part   #########################
##################################################################
hom=`pwd`
echo
echo "Welcome!"
echo "This is a script to judge the smae POSCAR files provided in the selected folder."
echo
if [[ "$1" != "" ]];then
  temp=$1
else
  echo
  echo "Current work directory :" ${hom}
  echo "==========================================="
  ls
  echo "==========================================="
  echo
  read -p "Which directory: " temp  #输入当前工作目录下的POSCAR文件所在目录
fi
if [[ ! -e "${hom}/${temp}" ]];then
  echo "No such file!"
  echo
  echo "Quit"
  echo "Have a good day!"
  exit 1
fi
work="${hom}/workout"
rm -rf ${work}
mkdir ${work}
hom="${hom}/${temp}"
pre=7
if [[ "$2" != "" ]];then
  pre=$2
else
  read -p "Precision: " temp
  if [ "${temp}" != "" ];then
    pre=${temp}
  fi
fi
echo
rm -rf ${hom}/compare_*
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
      phonopy -c PPOSCAR --symmetry>sym_${k}
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
    if [ "${judge[${order}]}" != "0" ];then
      continue
    fi
    
    # 比较一下POSCAR文件每个原子位置是不是完全一样的
    file1="POSCAR_g_${i}"
    file2="POSCAR_g_${j}"
    num_true=`awk 'END{print NR}' ${file1}`
    for ((k1=9;k1<=${num_true};k1++))
    do
      for l in 1 2 3
      do
        num1=`awk -v l=${l} 'NR=="'${k1}'"{print $l}' ${file1}`
        num2=`awk -v l=${l} 'NR=="'${k1}'"{print $l}' ${file2}`
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



改进版

```shell
#! /bin/bash
# posjudge
# This is a script to judge the smae POSCAR files provided in the selected folder.

##################################################################
##########################  Input part   #########################
##################################################################
hom=`pwd`
echo
echo "Welcome!"
echo "This is a script to judge the smae POSCAR files provided in the selected folder."
echo
if [[ "$1" != "" ]];then
  temp=$1
else
  echo
  echo "Current work directory :" ${hom}
  echo "==========================================="
  ls
  echo "==========================================="
  echo
  read -p "Which directory: " temp  #输入当前工作目录下的POSCAR文件所在目录
fi
if [[ ! -e "${hom}/${temp}" ]];then
  echo "No such file!"
  echo
  echo "Quit"
  echo "Have a good day!"
  exit 1
fi
work="${hom}/workout"
rm -rf ${work}
mkdir ${work}
hom="${hom}/${temp}"
pre=7
if [[ "$2" != "" ]];then
  pre=$2
else
  read -p "Precision: " temp
  if [ "${temp}" != "" ];then
    pre=${temp}
  fi
fi
echo
rm -rf ${hom}/compare_*
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
      phonopy -c PPOSCAR --symmetry>sym_${k} #注意：此时的sym文件是基于PPOSCAR的
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
    
    #晶体对称性比较
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
          for ((k2=k1+1;k2<=type_num[${type}];k2++))
          do
            eval w1='`'printf "%d" '\'"'""$"{atom${type}_Wyckoff_mark[$k1]}'`'
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
              elif (( x2 == x1 ));then
                eval y1='$'{atom${type}_irre[$k1]}
                eval y2='$'{atom${type}_irre[$k2]}
                if (( y2 < y1 ));then
                  eval temp='$'{atom${type}_irre['$'k1]}
                  eval atom${type}_irre['$'k1]='$'{atom${type}_irre['$'k2]}
                  eval atom${type}_irre['$'k2]='$'{temp}            
                  eval temp='$'{atom${type}_equiv_num['$'k1]}
                  eval atom${type}_equiv_num['$'k1]='$'{atom${type}_equiv_num['$'k2]}
                  eval atom${type}_equiv_num['$'k2]='$'{temp}       
                  eval temp='$'{atom${type}_Wyckoff_mark['$'k1]}
                  eval atom${type}_Wyckoff_mark['$'k1]='$'{atom${type}_Wyckoff_mark['$'k2]}
                  eval atom${type}_Wyckoff_mark['$'k2]='$'{temp}
                fi
              fi
            fi
          done
        done
      done
      # 生成用于比较的对称性信息文件sym_log_*,以及更加详细的版本sym_log_d_*文件（多一列个不等价原子在生成的PPOSCAR文件中的位置）
      touch ${dir}/sym_log_${ii}
      touch ${dir}/sym_log_d_${ii}
      echo ${group_num}>>${dir}/sym_log_${ii}
      echo ${group_num}>>${dir}/sym_log_d_${ii}
      for ((type=1;type<=${#atom_num[*]};type++))
      do
        for ((k=1;k<=type_num[${type}];k++))
        do      
          eval printf "%3s%4s%3s%3d'\n'" '$'{atom${type}_order['$'k]} '$'{atom_type['$'{type}]} '$'{atom${type}_Wyckoff_mark['$'k]} '$'{atom${type}_equiv_num['$'k]}|tee -a ${dir}/sym_log_${ii}
          eval printf "%3s%4s%3s%3d%3d'\n'" '$'{atom${type}_order['$'k]} '$'{atom_type['$'{type}]} '$'{atom${type}_Wyckoff_mark['$'k]} '$'{atom${type}_equiv_num['$'k]} '$'{atom${type}_irre['$'k]|tee -a ${dir}/sym_log_d_${ii}
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
    if [ "$judge[${order}]" != "0" ];then
      continue
    fi
    # 以上比较了PPOSCAR的晶格常数、原子总数、对称性信息（包含了原子种类及其个数）。如果说都是一样的话，那么接下来就要具体比较POSCAR文件了。策略：把POSCAR文件中首个元素中字母最小、重数最低的那类原子（对应于sym_log文件中的第一个原子）放到原点位置，比较两个POSCAR文件。如果遇到不等价原子数目一致的情况，按照原先POSCAR文件中的元素排列顺序来（要求两个POSCAR文件中元素排列顺序一致），如果遇到字母、重数都一样的情况，那么就固定其中的一个POSCAR，让另一个POSCAR文件把这些看起来信息一模一样的不等价位置都位移至原点尝试一下。
    
    # 确定要放到原点的原子对应的序号
    atom_origin=`awk 'NR==2{print $5}' ${dir}/sym_log_d_$i`
    # 读取原点原子位置信息
    x_move=
    y_move
    z_move
    
    
    
    
    
    
    
    
    
    
    
    
    
    # 比较一下POSCAR文件每个原子位置是不是完全一样的
    file1="POSCAR_g_${i}"
    file2="POSCAR_g_${j}"
    num_true=`awk 'END{print NR}' ${file1}`
    for ((k1=9;k1<=${num_true};k1++))
    do
      for l in 1 2 3
      do
        num1=`awk -v l=${l} 'NR=="'${k1}'"{print $l}' ${file1}`
        num2=`awk -v l=${l} 'NR=="'${k1}'"{print $l}' ${file2}`
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

