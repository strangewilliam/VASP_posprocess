# 整理POSCAR文件

## 1 初始化：除去Selective Dynamics 以及T/F标记，除去分子动力学的初速度

posinit

```shell
#! /bin/bash
# posinit
# 初始化：除去Selective Dynamics 以及T/F标记，除去分子动力学的初速度

##################################################################
##########################  Input part   #########################
##################################################################
hom=`pwd`
echo
echo "Welcome!"
echo "This is a script to initialize the POSCAR files provided."
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
  read -p "Which file: " temp  # 输入POSCAR文件名称
fi
if [[ ! -e "${hom}/${temp}" ]];then
  echo "No such file!"
  echo
  echo "Quit"
  echo "Have a good day!"
  exit 1
fi
file_name="${temp}" # POSCAR文件名称

if [[ "$2" != "" ]];then
  temp=$2
else
  read -p "New POSCAR file name: " temp #新的POSCAR文件命名
fi
if [ "${temp}" == "" ];then
  temp="${file_name}_i"
fi
output=${temp}
if [ -e "${hom}/${output}" ];then
  echo "Repeated file name!"
  echo "Finding appropriate NEW POSCAR file name ..."
  echo
  k=0
  while [ -e "${hom}/${output}_${k}" ]
  do
   let k++
  done
  output=${output}_${k}
fi
echo "The NEW POSCAR file name: ${output}"
echo
file_name="${hom}/${file_name}"
output=${hom}/${output}
touch ${output}  # 生成输出的POSCAR文件
atom_num=`awk 'NR==7{a=0;for(i=1;i<=NF;i++){a=a+$i};print a}' ${file_name}` # 原子总数
begin=`awk 'NR>7{j=toupper(substr($1,1,1));if(j=="D"||j=="C"){print NR}}' ${file_name}` #开始复制的行号
let end=begin+atom_num  # 结束复制的行号

##################################################################
##################################################################

##################################################################
##########################   Output part   #######################
##################################################################
# 前7行直接复制
head -7 ${file_name}>>${output}
# 后面的行只复制原子坐标
for i in $(seq ${begin} ${end})
do
  if ((i==begin));then
    awk -v j=${i} 'NR==j{print $0}' ${file_name} >> ${output}
  else
    awk -v j=${i} '{if(NR==j){printf ("%25.18f%25.18f%25.18f\n", $1,$2,$3)}}' ${file_name} >> ${output}
  fi
done
echo "Successfully generated NEW POSCAR file!"
echo "Finished."
echo "Have a good day!"
##################################################################
##################################################################
exit 0
```



## 2 调整POSCAR文件中元素整体的位置（也就是第六七两行中元素的位置）

​    个性化：自己输入确定排列位置【s】或者按照元素周期表的顺序排列【p】

posarrange

```shell
#! /bin/bash
# posarrange
# 调整POSCAR文件中元素整体的位置（也就是第六七两行中元素的位置）
# 个性化：自己输入确定排列位置【s】或者按照元素周期表的顺序排列【p】
# 必须要预先经过posinit处理的POSCAR文件才可以用这个脚本

##################################################################
##########################  Input part   #########################
##################################################################
hom=`pwd`
program="/fs12/home/zhj_xujq/bin" # 程序存放的默认位置
echo
echo "Welcome!"
echo "This is a script to rearrange the POSCAR files provided."
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
  read -p "Which file: " temp  # 输入POSCAR文件名称
fi
if [[ ! -e "${hom}/${temp}" ]];then
  echo "No such file!"
  echo
  echo "Quit"
  echo "Have a good day!"
  exit 1
fi
file_name="${temp}" # POSCAR文件名称

if [[ "$2" != "" ]];then
  temp=$2
else
  read -p "New POSCAR file name: " temp #新的POSCAR文件命名
fi
if [ "${temp}" == "" ];then
  temp="${file_name}_r"
fi
output=${temp}
if [ -e "${hom}/${output}" ];then
  echo "Repeated file name!"
  echo "Finding appropriate NEW POSCAR file name ..."
  echo
  k=0
  while [ -e "${hom}/${output}_${k}" ]
  do
    let k++
  done
  output=${output}_${k}
fi
echo "The NEW POSCAR file name: ${output}"
echo
file_name="${hom}/${file_name}"
output=${hom}/${output}
touch ${output}  # 生成输出的POSCAR文件

# 选择排序的方式：自己输入确定排列位置【s】或者按照元素周期表的顺序排列【p】
if [[ "$3" != "" ]];then
  temp=$3
else
  read -p "Arrange mode (Self[s] or Periodic table[p]):  " temp 
fi
if [[ "${temp}" != "s" && "${temp}" != "p" ]];then
  echo 
  echo "Wrong input!"
  echo "Quit."
  exit 1
fi
arrange=${temp} # 排序方式的标记

# 元素周期表的位置
if [[ "${arrange}" == "p" ]];then
  if [[ "$4" != "" ]];then
    program=$4
  else
    read -p "Where is the PERIODIC_TABLE? " temp
    if [[ "${temp}" != "" ]];then
      program=${temp}
    fi
  fi
fi

k=0
for temp in `awk 'NR==6{print $0}' ${file_name}`
do
  let k++
  atom_type[$k]=${temp}  # 原子种类
done
k=0
for temp in `awk 'NR==7{print $0}' ${file_name}`
do
  let k++
  type_num[$k]=${temp}  # 每种原子的个数
done

# 判断*POSCAR*文件是否符合要求
judge=`awk 'NR==8{j=toupper(substr($1,1,1));if(j=="D"||j=="C"){print NR}}' ${file_name}`
if [[ "${judge}" != "8" ]];then
  echo "Wrong POSCAR file format!"
  echo "Quit."
  rm ${output}
  exit 1
fi

begin=9 #开始复制的行号
for ((i=1;i<=${#atom_type[*]};i++))
do
  atom_begin[$i]=${begin}
  let atom_end[$i]=atom_begin[$i]+type_num[$i]-1 # 每个元素起始行号
  let begin=atom_end[$i]+1 # 每个元素终止行号
done

##################################################################
##################################################################

##################################################################
######################   Data process part   #####################
##################################################################
# 生成排序数组order[*]
case ${arrange} in
  "s")
    # 生成ARRANGE文件。将来可以对ARRANGE文件的格式进行说明
    while [ ! -e "${hom}/ARRANGE" ]
    do
      echo
      echo "POSCAR file: "
      echo "==========================================="
      head -7 ${file_name}|tail -2
      echo "==========================================="
      echo
      temp=""
      while [ "${temp}" == "" ]
      do
        read -p "The order of elements: " temp
      done
      touch ${hom}/ARRANGE
      for things in ${temp}
      do
        echo ${things} >>${hom}/ARRANGE
      done
    done
    i=0
    for temp1 in ${atom_type[*]}
    do
      let i++
      j=1
      for temp2 in `awk '{print $1}' ${hom}/ARRANGE`
      do
        if [ "${temp1}" == "${temp2}" ];then
          order[$i]=$j  # 生成元素排序
          order_num[$i]=$i  # 元素对应序号
          break
        fi
        let j++
      done
      if ((j>${#atom_type[*]}));then # 判断ARRANGE文件是否正确
        echo "Wrong ARRANGE file!"
        echo "Quit."
        rm ${output}
        exit 1
      fi
    done
  ;;
  "p")
    if [ ! -e "${program}/PERIODIC_TABLE" ];then
      echo "*PERIODIC_TABLE* file not found. Please provide it!"
      echo "Quit."
      rm ${output}
      exit 1
    fi
    i=0
    for temp1 in ${atom_type[*]}
    do
      let i++
      j=1
      for temp2 in `awk '{print $1}' ${program}/PERIODIC_TABLE`
      do
        if [ "${temp1}" == "${temp2}" ];then
          order[$i]=$j  # 生成元素排序 (元素周期表中的序号)
          order_num[$i]=$i  # 元素对应序号
          break
        fi
        let j++
      done
      if ((j>118));then # 判断POSCAR文件是否正确
        echo "Wrong element *${temp1}* provided in POSCAR file!"
        echo "Quit."
        rm ${output}
        exit 1
      fi
    done
  ;;
esac

# 排序
for ((k1=1;k1<=${#order[*]}-1;k1++))
do 
  for ((k2=k1+1;k2<=${#order[*]};k2++))
  do
    if (( ${order[${k1}]} > ${order[${k2}]} ));then
      temp=${order[${k1}]}
      order[${k1}]=${order[${k2}]}
      order[${k2}]=${temp}
      temp=${order_num[${k1}]}
      order_num[${k1}]=${order_num[${k2}]}
      order_num[${k2}]=${temp}      
    fi
  done
done

##################################################################
##################################################################

##################################################################
##########################   Output part   #######################
##################################################################
# 前5行直接复制
head -5 ${file_name}>>${output}
# 后面的行要先排序后复制
# 修改第6、7行
format="";string6="";string7=""
for k in ${order_num[*]}
do
  format="${format}%5s"
  string6="${string6} ${atom_type[$k]}"
  string7="${string7} ${type_num[$k]}"
done  
eval printf '"'${format}'\'n'"' ${string6}>>${output}
eval printf '"'${format}'\'n'"' ${string7}>>${output}
# 第8行直接复制
head -8 ${file_name}|tail -1>>${output}
for k in ${order_num[*]}
do
  for i in $(seq ${atom_begin[$k]} ${atom_end[$k]})
  do
    awk -v j=${i} '{if(NR==j){printf ("%25.18f%25.18f%25.18f\n", $1,$2,$3)}}' ${file_name} >> ${output}
  done
done
echo "Successfully generated NEW POSCAR file!"
echo "Finished."
echo "Have a good day!"
##################################################################
##################################################################
exit 0
```

## 3 把POSCAR文件中一些负值或者比较接近1的坐标移回到原点0位置（需要设置一定精度）

posadjust

```shell
#! /bin/bash
# posadjust
# 把POSCAR文件中一些负值或者比较接近或者超过1的坐标移回到原点0位置（需要设置一定精度）
# 要求POSCAR文件是Direct坐标
# 需要设置精度 此程序默认精度为1E-7
# 不修改Selective Dynamics

##################################################################
##########################  Input part   #########################
##################################################################
hom=`pwd`
echo
echo "Welcome!"
echo "This is a script to transfer some coordinates that are close to or exceed 1 back to the range 0~1."
echo "The POSCAR file format must be *Direct coordinate*."
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
  read -p "Which file: " temp  # 输入POSCAR文件名称
fi
if [[ ! -e "${hom}/${temp}" ]];then
  echo "No such file!"
  echo
  echo "Quit"
  echo "Have a good day!"
  exit 1
fi
file_name="${temp}" # POSCAR文件名称

if [[ "$2" != "" ]];then
  temp=$2
else
  read -p "New POSCAR file name: " temp #新的POSCAR文件命名
fi
if [ "${temp}" == "" ];then
  temp="${file_name}_a"
fi
output=${temp}
if [ -e "${hom}/${output}" ];then
  echo "Repeated file name!"
  echo "Finding appropriate NEW POSCAR file name ..."
  echo
  k=0
  while [ -e "${hom}/${output}_${k}" ]
  do
    let k++
  done
  output=${output}_${k}
fi
echo "The NEW POSCAR file name: ${output}"
echo
file_name="${hom}/${file_name}"
output=${hom}/${output}
touch ${output}  # 生成输出的POSCAR文件

# 判断数据相等的精度 
if [[ "$3" != "" ]];then
  prec=$3
else
  read -p "Precision[Default= 1E-7]: " temp
  if [ "${temp}" != "" ];then
    prec=${temp}
  else
    prec="1E-7" # 默认精度为1E-7
  fi
fi
echo "Precision: ${prec}"
echo
p=${prec:0-1:1}
let pp=p+2
prec=`printf "%${pp}.${p}f\n" ${prec}`
k=0
for temp in `awk 'NR==7{print $0}' ${file_name}`
do
  let k++
  type_num[$k]=${temp}  # 每种原子的个数
done

atom_num=`awk 'NR==7{a=0;for(i=1;i<=NF;i++){a=a+$i};print a}' ${file_name}` # 原子总数
begin=`awk 'NR>7{j=toupper(substr($1,1,1));if(j=="D"||j=="C"){print NR+1}}' ${file_name}` #开始处理的行号
let end=begin+atom_num-1  # 结束处理的行号
##################################################################
##################################################################

##################################################################
##########################   Output part   #######################
##################################################################
# 前*begin-1*行直接复制
let temp=begin-1
head -${temp} ${file_name}>>${output}
# 从*begin*行开始后面的行进行数据处理
# 首先对该坐标对1取余数
# 接着判断该坐标是否和1在精度prec范围内比较接近，如果是，那么将其改成0
# 逐行处理
for k in $(seq ${begin} ${end})
do
  temp1=`awk 'NR=="'$k'"{print $1}' ${file_name}`
  temp1=`echo "${temp1}%1"|bc`
  temp=${temp1:0:1}
  if [ ${temp} == "-" ];then
    temp1=`echo "(1+${temp1})"|bc`
  fi
  judge=`echo "(1-${temp1})/${prec}"|bc`
  if ((judge<=1));then
    temp1=0
  fi
  temp2=`awk 'NR=="'$k'"{print $2}' ${file_name}`
  temp2=`echo "${temp2}%1"|bc`
  temp=${temp2:0:1}
  if [ ${temp} == "-" ];then
    temp2=`echo "(1+${temp2})"|bc`
  fi
  judge=`echo "(1-${temp2})/${prec}"|bc`
  if ((judge<=1));then    
    temp2=0
  fi
  temp3=`awk 'NR=="'$k'"{print $3}' ${file_name}`
  temp3=`echo "${temp3}%1"|bc`
  temp=${temp3:0:1}
  if [ ${temp} == "-" ];then
    temp3=`echo "(1+${temp3})"|bc`
  fi
  judge=`echo "(1-${temp3})/${prec}"|bc`
  if ((judge<=1));then 
    temp3=0
  fi
  temp=`awk 'NR==8{j=toupper(substr($1,1,1));print (j=="S")?(1):(-1)}' ${file_name}` #判断是否为Selective Dynamics 模式
  if [[ "${temp}" == "1" ]];then
    temp4=`awk 'NR=="'$k'"{print $4}' ${file_name}`
    temp5=`awk 'NR=="'$k'"{print $5}' ${file_name}`
    temp6=`awk 'NR=="'$k'"{print $6}' ${file_name}`
    printf "%25.18f%25.18f%25.18f%2s%2s%2s\n" ${temp1} ${temp2} ${temp3} ${temp4} ${temp5} ${temp6}>>${output}
  else
    printf "%25.18f%25.18f%25.18f\n" ${temp1} ${temp2} ${temp3}>>${output}
  fi
done

echo "Successfully generated NEW POSCAR file!"
echo "Finished."
echo "Have a good day!"
##################################################################
##################################################################
exit 0
```

## 4 每个元素坐标按照x, y, z 顺序依次升序排列

posorder

```shell
#! /bin/bash
# posorder
# 把POSCAR文件中每个元素坐标按照x, y, z 顺序依次升序排列
# 不修改Selective Dynamics
# 需要提供相等的精度prec

##################################################################
########################  function part   ########################
##################################################################
function same(){
  ppre=$3
  ppre=${ppre:0-1:1}
  pre=1
  for ((fj=1;fj<=${ppre};fj++))
  do
    pre="${pre}0"
  done   
  abs=`awk -v p=${pre} -v n1=$1 -v n2=$2 'BEGIN{a=p*(n1-n2);print a}'` 
  abs=`awk -v abs=${abs} 'BEGIN{print (abs<0)?(-abs):(abs)}'`
  judge_t=`awk -v abs=${abs} 'BEGIN{print (abs<1)?0:1}'`
  return ${judge_t}
}
##################################################################
##################################################################


##################################################################
##########################  Input part   #########################
##################################################################
hom=`pwd`
echo
echo "Welcome!"
echo "This is a script to arrange the coordinate in order"
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
  read -p "Which file: " temp  # 输入POSCAR文件名称
fi
if [[ ! -e "${hom}/${temp}" ]];then
  echo "No such file!"
  echo
  echo "Quit"
  echo "Have a good day!"
  exit 1
fi
file_name="${temp}" # POSCAR文件名称

if [[ "$2" != "" ]];then
  temp=$2
else
  read -p "New POSCAR file name: " temp #新的POSCAR文件命名
fi
if [ "${temp}" == "" ];then
  temp="${file_name}_o"
fi
output=${temp}
if [ -e "${hom}/${output}" ];then
  echo "Repeated file name!"
  echo "Finding appropriate NEW POSCAR file name ..."
  echo
  k=0
  while [ -e "${hom}/${output}_${k}" ]
  do
    let k++
  done
  output=${output}_${k}
fi
echo "The NEW POSCAR file name: ${output}"
echo
file_name="${hom}/${file_name}"
output=${hom}/${output}
touch ${output}  # 生成输出的POSCAR文件

# 判断数据相等的精度 
if [[ "$3" != "" ]];then
  prec=$3
else
  read -p "Precision[Default= 1E-7]: " temp
  if [ "${temp}" != "" ];then
    prec=${temp}
  else
    prec="1E-7" # 默认精度为1E-7
  fi
fi
echo "Precision: ${prec}"
echo

k=0
for temp in `awk 'NR==7{print $0}' ${file_name}`
do
  let k++
  type_num[$k]=${temp}  # 每种原子的个数
done

atom_num=`awk 'NR==7{a=0;for(i=1;i<=NF;i++){a=a+$i};print a}' ${file_name}` # 原子总数
begin=`awk 'NR>7{j=toupper(substr($1,1,1));if(j=="D"||j=="C"){print NR+1}}' ${file_name}` #开始处理的行号
let end=begin+atom_num-1  # 结束处理的行号

temp=${begin}
for ((i=1;i<=${#type_num[*]};i++))
do
  atom_begin[$i]=${temp}
  let atom_end[$i]=atom_begin[$i]+type_num[$i]-1 # 每个元素起始行号
  let temp=atom_end[$i]+1 # 每个元素终止行号
done

k=0  # order_num[*]数组中的位置
for i in $(seq ${begin} ${end})
do
  let k++
  order_num[$k]=$i  # 每个原子对应行号
done
##################################################################
##################################################################

##################################################################
##########################   Output part   #######################
##################################################################
# 前*begin-1*行直接复制
let temp=begin-1
head -${temp} ${file_name}>>${output}
# 从*begin*行开始后面的行进行数据处理
# 写入第type种元素原子坐标到temp*
for ((i=1;i<=atom_num;i++))
do
  temp1[$i]=`awk -v l=${order_num[$i]} 'NR==l{print $1}' ${file_name}`
  temp2[$i]=`awk -v l=${order_num[$i]} 'NR==l{print $2}' ${file_name}`
  temp3[$i]=`awk -v l=${order_num[$i]} 'NR==l{print $3}' ${file_name}`
done
# 排序
add=1
for ((type=1;type<=${#type_num[*]};type++))
do
  for ((i=add;i<=add+type_num[${type}]-2;i++))
  do
    for ((j=i+1;j<=add+type_num[${type}]-1;j++))
    do
      # 判断谁大谁小
      #先判断是否在一定精度范围内相等
      same ${temp1[$i]} ${temp1[$j]} ${prec}
      judge_same=$?
      judge=`awk -v t1=${temp1[$i]} -v t2=${temp1[$j]} 'BEGIN{d=t1-t2;print (d>0)?(1):(-1)}'`
      if ((judge==1&&judge_same==1));then
        temp=${temp1[$i]}
        temp1[$i]=${temp1[$j]}
        temp1[$j]=${temp}
        temp=${temp2[$i]}
        temp2[$i]=${temp2[$j]}
        temp2[$j]=${temp}
        temp=${temp3[$i]}
        temp3[$i]=${temp3[$j]}
        temp3[$j]=${temp}
        temp=${order_num[$i]}
        order_num[$i]=${order_num[$j]}
        order_num[$j]=${temp}
      elif ((judge_same==0));then
        #先判断是否在一定精度范围内相等
        same ${temp2[$i]} ${temp2[$j]} ${prec}
        judge_same=$?
        judge=`awk -v t1=${temp2[$i]} -v t2=${temp2[$j]} 'BEGIN{d=t1-t2;print (d>0)?(1):(-1)}'`
        if ((judge==1&&judge_same==1));then
          temp=${temp1[$i]}
          temp1[$i]=${temp1[$j]}
          temp1[$j]=${temp}
          temp=${temp2[$i]}
          temp2[$i]=${temp2[$j]}
          temp2[$j]=${temp}
          temp=${temp3[$i]}
          temp3[$i]=${temp3[$j]}
          temp3[$j]=${temp}
          temp=${order_num[$i]}
          order_num[$i]=${order_num[$j]}
          order_num[$j]=${temp}
        elif ((judge_same==0));then
          #先判断是否在一定精度范围内相等
          same ${temp3[$i]} ${temp3[$j]} ${prec}
          judge_same=$?
          judge=`awk -v t1=${temp3[$i]} -v t2=${temp3[$j]} 'BEGIN{d=t1-t2;print (d>0)?(1):(-1)}'`
          if ((judge==1&&judge_same==1));then
            temp=${temp1[$i]}
            temp1[$i]=${temp1[$j]}
            temp1[$j]=${temp}
            temp=${temp2[$i]}
            temp2[$i]=${temp2[$j]}
            temp2[$j]=${temp}
            temp=${temp3[$i]}
            temp3[$i]=${temp3[$j]}
            temp3[$j]=${temp}
            temp=${order_num[$i]}
            order_num[$i]=${order_num[$j]}
            order_num[$j]=${temp}
          elif ((judge_same==0));then
            echo "The same two atom ${i} and ${j} appear in the POSCAR file! Wrong!"
            echo "Quit."
            rm ${output}
            exit 1
          fi
        fi
      fi
    done
  done
  let add=add+type_num[${type}]
done


for ((k=1;k<=atom_num;k++))
do
  printf "%25.18f%25.18f%25.18f\n" ${temp1[$k]} ${temp2[$k]} ${temp3[$k]}>>${output}
done


#for k in ${order_num[*]}
#do
#  awk -v j=${k} 'NR==j{print $0}' ${file_name} >> ${output}
#done

echo "Successfully generated NEW POSCAR file!"
echo "Finished."
echo "Have a good day!"
##################################################################
##################################################################
exit 0
```

## 5 给定一个平移矢量，改动POSCAR文件中所有原子的位置信息

postrans

需要设置精度 此程序默认精度为1E-7

```shell
#! /bin/bash
# postrans
# 给定一个平移矢量，改动POSCAR文件中所有原子的位置信息,挪动的是原点的位置。所以新的坐标是减去平移矢量
# 需要设置精度 此程序默认精度为1E-7
# 不修改Selective Dynamics
# 该脚本拥有*posadjust*的所有功能，可以完全替代
# postrans POSCAR POSCAR_t 1E-7 0.5 0.5 0.5
##################################################################
##########################  Input part   #########################
##################################################################
hom=`pwd`
echo
echo "Welcome!"
echo "1. This is a script to change all coordinates of all atoms with respect to the translation vector given."
echo "2. The POSCAR file format must be *Direct coordinate*."
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
  read -p "Which file: " temp  # 输入POSCAR文件名称
fi
if [[ ! -e "${hom}/${temp}" ]];then
  echo "No such file!"
  echo
  echo "Quit"
  echo "Have a good day!"
  exit 1
fi
file_name="${temp}" # POSCAR文件名称

if [[ "$2" != "" ]];then
  temp=$2
else
  read -p "New POSCAR file name: " temp #新的POSCAR文件命名
fi
if [ "${temp}" == "" ];then
  temp="${file_name}_t"
fi
output=${temp}
if [ -e "${hom}/${output}" ];then
  echo "Repeated file name!"
  echo "Finding appropriate NEW POSCAR file name ..."
  echo
  k=0
  while [ -e "${hom}/${output}_${k}" ]
  do
    let k++
  done
  output=${output}_${k}
fi
echo "The NEW POSCAR file name: ${output}"
echo
file_name="${hom}/${file_name}"
output=${hom}/${output}
touch ${output}  # 生成输出的POSCAR文件

# 判断数据相等的精度 
if [[ "$3" != "" ]];then
  prec=$3
else
  read -p "Precision[Default= 1E-7]: " temp
  if [ "${temp}" != "" ];then
    prec=${temp}
  else
    prec="1E-7" # 默认精度为1E-7
  fi
fi
echo "Precision: ${prec}"
echo
p=${prec:0-1:1}
let pp=p+2
prec=`printf "%${pp}.${p}f\n" ${prec}`

# 输入平移矢量的坐标
if [[ "$4" == "" || "$5" == "" || "$6" == "" ]];then
  read -p "Translation vevtor[Direct coordinate format]: " temp
  if [ "temp" == "" ];then
    temp="0 0 0"
  fi
else
  temp="$4 $5 $6"
fi
k=0
for i in ${temp}
do
  let k++
  translate[$k]=$i  # 平移矢量的坐标存放
done
echo "translation vector:"
printf "%10.7f\n%10.7f\n%10.7f\n" ${translate[1]} ${translate[2]} ${translate[3]}
if [ "${#translate[*]}" != "3" ];then
  echo "Wrong translation vector input!"
  echo "Quit."
  rm ${output}
  exit 1
fi

k=0
for temp in `awk 'NR==7{print $0}' ${file_name}`
do
  let k++
  type_num[$k]=${temp}  # 每种原子的个数
done

atom_num=`awk 'NR==7{a=0;for(i=1;i<=NF;i++){a=a+$i};print a}' ${file_name}` # 原子总数
begin=`awk 'NR>7{j=toupper(substr($1,1,1));if(j=="D"||j=="C"){print NR+1}}' ${file_name}` #开始处理的行号
let end=begin+atom_num-1  # 结束处理的行号
##################################################################
##################################################################

##################################################################
##########################   Output part   #######################
##################################################################
# 前*begin-1*行直接复制
let temp=begin-1
head -${temp} ${file_name}>>${output}
# 从*begin*行开始后面的行进行数据处理
# 首先进行坐标平移，逐行处理
# 接着对该坐标对1取余数
# 接着判断该坐标是否和1在精度prec范围内比较接近，如果是，那么将其改成0

for k in $(seq ${begin} ${end})
do
  temp1=`awk 'NR=="'$k'"{print $1}' ${file_name}`
  temp1=`echo "(${temp1}-(${translate[1]}))%1"|bc` # 添加平移，接着取余数
  temp=${temp1:0:1}
  if [ "${temp}" == "-" ];then
    temp1=`echo "(1+${temp1})"|bc`
  fi
  judge=`echo "(1-${temp1})/${prec}"|bc`
  if ((judge<=1));then
    temp1=0
  fi
  temp2=`awk 'NR=="'$k'"{print $2}' ${file_name}`
  temp2=`echo "(${temp2}-(${translate[2]}))%1"|bc`  # 添加平移，接着取余数
  temp=${temp2:0:1}
  if [ "${temp}" == "-" ];then
    temp2=`echo "(1+${temp2})"|bc`
  fi
  judge=`echo "(1-${temp2})/${prec}"|bc`
  if ((judge<=1));then    
    temp2=0
  fi
  temp3=`awk 'NR=="'$k'"{print $3}' ${file_name}`
  temp3=`echo "(${temp3}-(${translate[3]}))%1"|bc`  # 添加平移，接着取余数
  temp=${temp3:0:1}
  if [ "${temp}" == "-" ];then
    temp3=`echo "(1+${temp3})"|bc`
  fi
  judge=`echo "(1-${temp3})/${prec}"|bc`
  if ((judge<=1));then 
    temp3=0
  fi
  temp=`awk 'NR==8{j=toupper(substr($1,1,1));print (j=="S")?(1):(-1)}' ${file_name}` #判断是否为Selective Dynamics 模式
  if [[ "${temp}" == "1" ]];then
    temp4=`awk 'NR=="'$k'"{print $4}' ${file_name}`
    temp5=`awk 'NR=="'$k'"{print $5}' ${file_name}`
    temp6=`awk 'NR=="'$k'"{print $6}' ${file_name}`
    printf "%25.18f%25.18f%25.18f%2s%2s%2s\n" ${temp1} ${temp2} ${temp3} ${temp4} ${temp5} ${temp6}>>${output}
  else
    printf "%25.18f%25.18f%25.18f\n" ${temp1} ${temp2} ${temp3}>>${output}
  fi
done

echo "Successfully generated NEW POSCAR file!"
echo "Finished."
echo "Have a good day!"
##################################################################
##################################################################
exit 0
```



## 6 给定一个方向矢量以及旋转角度，改动POSCAR文件中所有原子的位置信息（大业未成）

```shell
#! /bin/bash

##################################################################
########################  constant part   ########################
##################################################################
pi=3.1415926535897932384626433
function same(){
  ppre=$3
  ppre=${ppre:0-1:1}
  pre=1
  for ((fj=1;fj<=${ppre};fj++))
  do
    pre="${pre}0"
  done   
  abs=`awk -v p=${pre} -v n1=$1 -v n2=$2 'BEGIN{a=p*(n1-n2);print a}'` 
  abs=`awk -v abs=${abs} 'BEGIN{print (abs<0)?(-abs):(abs)}'`
  judge_t=`awk -v abs=${abs} 'BEGIN{print (abs<1)?0:1}'`
  return ${judge_t}
}
##################################################################
##################################################################

```



## 7 计算给定原子与所有其他原子的距离

posdist

```shell
#! /bin/bash
# posdist
# 选定一个原子，计算其与其他所有晶胞内原子之间的距离，并按照不同元素排序输出至文档
# 需要前面的postrans脚本的辅助
# 需要设置精度 此程序默认精度为1E-7,最高精度为1E-9,最低精度为1E0
# posdist POSCAR 18 POSCAR_dist 1E-7
##################################################################
########################  Function part   ########################
##################################################################
function distance(){
  local temp=`echo "scale=$7;sqrt(($1-($4))^2+($2-($5))^2+($3-($6))^2)"|bc`

  printf "%25.18f\n" ${temp}
  return 0
}

function same(){
  local ppre=$3
  ppre=${ppre:0-1:1}
  local pre=1
  for ((fj=1;fj<=${ppre};fj++))
  do
    pre="${pre}0"
  done   
  local abs=`awk -v p=${pre} -v n1=$1 -v n2=$2 'BEGIN{a=p*(n1-n2);print a}'` 
  abs=`awk -v abs=${abs} 'BEGIN{print (abs<0)?(-abs):(abs)}'`
  local judge_t=`awk -v abs=${abs} 'BEGIN{print (abs<1)?0:1}'`
  return ${judge_t}
}
##################################################################
##################################################################

##################################################################
##########################  Input part   #########################
##################################################################
hom=`pwd`
echo
echo "*posdist* script"
echo
if [[ "$1" != "" ]];then
  temp=$1
else
  echo "Welcome!"
  echo "1. This is a script to calculate the distance between one selected atom and all of the other atoms in the crystal cell."
  echo "2. The script *postrans* are needed."
  echo
  echo "Current work directory :" ${hom}
  echo "==========================================="
  ls
  echo "==========================================="
  echo
  read -p "Which file: " temp  # 输入POSCAR文件名称
fi
if [[ ! -e "${hom}/${temp}" ]];then
  echo "No such file!"
  echo
  echo "Quit"
  echo "Have a good day!"
  exit 1
fi
file_name="${temp}" # POSCAR文件名称

if [[ "$2" != "" ]];then
  temp=$2
else
  echo "POSCAR file: "
  echo "==========================================="
  head -7 ${file_name}|tail -2
  echo "==========================================="
  echo
  temp=""
  while [ "${temp}" == "" ]
  do
    read -p "Which atom[e.g. 4]: " temp
  done
fi
atom_selected=${temp}  # 选定要计算距离的那一个原子对应的序号

scale=`awk 'NR==2{print $1}' ${file_name}` # The universal scaling factor (lattice constant), which is used to scale all lattice vectors and all atomic coordinates (of this value is negative it is interpreted as the total volume of the cell).
# 判断POSCAR文件是否符合要求
judge=`awk -v s=${scale} 'BEGIN{print (s>0)?(1):(-1)}'`
if ((judge<0));then
  echo "Wrong POSCAR format: the universal scaling factor on line 2 in POSCAR file should be POSITIVE."
  echo "No such file!"
  echo
  echo "Quit"
  echo "Have a good day!"
  exit 1
fi

begin=`awk 'NR>7{j=toupper(substr($1,1,1));if(j=="D"||j=="C"){print NR+1}}' ${file_name}` #开始处理的行号
let end=begin+atom_num-1  # 结束处理的行号
# 判断POSCAR格式是Direct还是Cartesian,并输入晶格矢量坐标
f_judge=`awk -v b=${begin} 'NR==b-1{j=toupper(substr($1,1,1));print j}' ${file_name}`
case ${f_judge} in
"D")
  echo "The POSCAR file is in *Direct* coordinate format."
  for k in 1 2 3
  do
    a_axis[$k]=`awk -v s=${scale} -v k=$k 'NR==3{a=s*$k;printf "%23.16f",a}' ${file_name}`
    b_axis[$k]=`awk -v s=${scale} -v k=$k 'NR==4{a=s*$k;printf "%23.16f",a}' ${file_name}`
    c_axis[$k]=`awk -v s=${scale} -v k=$k 'NR==5{a=s*$k;printf "%23.16f",a}' ${file_name}`
  done  
;;
"C")
  echo "The POSCAR file is in *Cartesian* coordinate format."
  for k in 1 2 3
  do
    a_axis[$k]=${scale};b_axis[$k]=${scale};c_axis[$k]=${scale}
  done
;;
esac

k=0
for temp in `awk 'NR==6{print $0}' ${file_name}`
do
  let k++
  atom_type[$k]=${temp}  # 原子种类
done

k=0;atom_num=0
for temp in `awk 'NR==7{print $0}' ${file_name}`
do
  let k++
  type_num[$k]=${temp}  # 每种原子的个数
  let total_temp=atom_num+type_num[$k]
  if ((atom_selected>atom_num && atom_selected<=total_temp));then
    let minus_temp=atom_selected-atom_num
    first_line="${atom_type[$k]}_${minus_temp}" # 输出文件第一行为选定元素的第几个原子
    echo "The selected atom: ${first_line}" # 输出选择的那一种元素的第几个原子
    echo
  fi
  atom_num=${total_temp} # 原子总数
done

let temp=atom_selected+begin-1
echo "The coordinates:" 
head -${temp} ${file_name}|tail -1 # 输出选择的那个原子的坐标信息
echo
temp=${begin}
for ((i=1;i<=${#type_num[*]};i++))
do
  atom_begin[$i]=${temp} # 每个元素起始行号
  let atom_end[$i]=atom_begin[$i]+type_num[$i]-1 # 每个元素终止行号
  let temp=atom_end[$i]+1
done


if [[ "$3" != "" ]];then
  temp=$3
else
  read -p "Output file name: " temp # 距离信息输出文件命名
fi
if [ "${temp}" == "" ];then
  temp="${file_name}_dist"
fi
output=${temp}
if [ -e "${hom}/${output}" ];then
  echo "Repeated file name!"
  echo "Finding appropriate NEW output file name ..."
  echo
  k=0
  while [ -e "${hom}/${output}_${k}" ]
  do
    let k++
  done
  output=${output}_${k}
fi
echo "The output file name: ${output}"
echo
echo ${first_line} > ${hom}/${output}  # 生成输出的POSCAR文件

# 判断数据相等的精度 
if [[ "$4" != "" ]];then
  prec=$4
else
  read -p "Precision[Default= 1E-7]: " temp
  if [ "${temp}" != "" ];then
    prec=${temp}
  else
    prec="1E-7" # 默认精度为1E-7
  fi
fi
echo "Precision: ${prec}"
echo
##################################################################
##################################################################

##################################################################
######################   Data process part   #####################
##################################################################
let temp=atom_selected+begin-1
# 输入选定原子的相对坐标信息
k=0
for things in `awk -v t=${temp} 'NR==t{print $1,$2,$3}' ${file_name}`
do
  let k++
  r_r[$k]=${things}  # 相对坐标
  t[$k]=`echo "${things}-0.5"|bc` # 把选定原子挪到晶胞的中间的平移矢量
done
# 把选定原子挪到晶胞的中间
rm -rf ${file_name}_t
postrans ${file_name} ${file_name}_t ${prec} ${t[1]} ${t[2]} ${t[3]}>/dev/null
file_name="${file_name}_t"  # 接下来就是处理 *_t 文件了
# 输入选定原子的笛卡尔坐标信息
for k in 1 2 3
do
  r_s[$k]=`echo "scale=20;0.5*${a_axis[$k]}+0.5*${b_axis[$k]}+0.5*${c_axis[$k]}"|bc` # 绝对坐标
done

# 计算距离并输出到文件${output}中
for ((j=1;j<=${#type_num[*]};j++))
do
  num_temp=0 # 记录某元素的第几个原子
  for ((i=atom_begin[$j];i<=atom_end[$j];i++))
  do
    let num_temp++
    k=0
    for things in `awk -v t=$i 'NR==t{print $1,$2,$3}' ${file_name}`
    do
      let k++
      r_r[$k]=${things}  #  每个原子的相对坐标信息
    done
    # 每个原子的笛卡尔坐标信息
    for k in 1 2 3
    do
      r_t[$k]=`echo "scale=20;${r_r[1]}*${a_axis[$k]}+${r_r[2]}*${b_axis[$k]}+${r_r[3]}*${c_axis[$k]}"|bc`
    done
    dist_temp=`distance ${r_s[1]} ${r_s[2]} ${r_s[3]} ${r_t[1]} ${r_t[2]} ${r_t[3]} 20`
    printf "%-8s%25.18f\n" ${atom_type[$j]}_${num_temp} ${dist_temp}>>${output}
  done
done
##################################################################
##################################################################

##################################################################
##########################   Output part   #######################
##################################################################
rm -rf ${hom}/${file_name}
echo "Successfully generated the distance file *${output}*!"
echo "Finished."
echo "Have a good day!"
##################################################################
##################################################################
exit 0
```



## 8 给定坐标变换矩阵以及平移矢量，改变坐标系（还需要进一步调试）

posrotrans

```shell
#! /bin/bash
# posrotrans
# 给定坐标变换矩阵以及平移矢量，改变坐标系
# 给定一个平移矢量，改动POSCAR文件中所有原子的位置信息,挪动的是原点的位置。所以新的坐标是减去平移矢量
# 需要设置精度 此程序默认精度为1E-7
# 不修改Selective Dynamics
# 该脚本拥有*posadjust*、*postrans*的所有功能，可以完全替代
# 需要提供rotr.in 文件，文件格式：
# rotr.in
# ===========
# 1/2 0 1
# 0 2 1
# 1 0.3 1
# 0.5 0.5 0.5
# ===========
# 注释：前三行为旋转矩阵的逆矩阵，最后一行为平移矢量
# 注意：旋转矩阵的逆矩阵行列式的绝对值必须为1
# 格式：posrotrans POSCAR POSCAR_rt 1E-7
# The rotation matrix is P, and the translation vector is T
# The old origin is (0,0,0)', and the new origin (q10,q20.q30)'=T
# The old coordinate is (r1,r2,r3)', and the new coordinate (p1,p2,p3)'=(r1,r2,r3)'-T
# The old basis vectors are (a1,a2,a3), and the new basis vectors are (b1,b2,b3)
# The basis vectors transformation (b1,b2,b3)=(a1,a2,a3)P
# (i,j,k)(x1,x2,x3)'=(a1,a2,a3)(p1,p2,p3)'
#                   =(b1,b2,b3)P^{-1}(p1,p2,p3)'
# (q1,q2,q3)'=P^{-1}(p1,p2,p3)'
#            =P^{-1}[(r1,r2,r3)'-T]
##################################################################
##########################  Input part   #########################
##################################################################
hom=`pwd`
echo
echo "Welcome!"
echo "1. This is a script to change all coordinates of all atoms with respect to the rotation matrix and the translation vector given in the *rotr.in* file."
echo "2. The POSCAR file format must be *Direct coordinate*."
echo
echo "The format of rotr.in"
echo "==========="
echo "1/2 0 1"
echo "0 2 1"
echo "1 0.3 1"
echo "0.5 0.5 0.5"
echo "==========="
echo 

# Judge if the rotr.in file exists
if  [ ! -e "rotr.in" ];then
  let judge++
  echo "${judge}. Please provide the rotr.in file!"
  echo "Have a good day!"
  exit 1
fi

if [[ "$1" != "" ]];then
  temp=$1
else
  echo
  echo "Current work directory :" ${hom}
  echo "==========================================="
  ls
  echo "==========================================="
  echo
  read -p "Which file: " temp  # 输入POSCAR文件名称
fi
if [[ ! -e "${hom}/${temp}" ]];then
  echo "No such file!"
  echo
  echo "Quit"
  echo "Have a good day!"
  exit 1
fi
file_name="${temp}" # POSCAR文件名称

if [[ "$2" != "" ]];then
  temp=$2
else
  read -p "New POSCAR file name: " temp #新的POSCAR文件命名
fi
if [ "${temp}" == "" ];then
  temp="${file_name}_rt"
fi
output=${temp}
if [ -e "${hom}/${output}" ];then
  echo "Repeated file name!"
  echo "Finding appropriate NEW POSCAR file name ..."
  echo
  k=0
  while [ -e "${hom}/${output}_${k}" ]
  do
    let k++
  done
  output=${output}_${k}
fi
echo "The NEW POSCAR file name: ${output}"
echo
#file_name="${hom}/${file_name}"
#output=${hom}/${output}
touch ${output}  # 生成输出的POSCAR文件

# 判断数据相等的精度 
if [[ "$3" != "" ]];then
  prec=$3
else
  read -p "Precision[Default= 1E-7]: " temp
  if [ "${temp}" != "" ];then
    prec=${temp}
  else
    prec="1E-7" # 默认精度为1E-7
  fi
fi
echo "Precision: ${prec}"
echo
p=${prec:0-1:1}
let pp=p+2
prec=`printf "%${pp}.${p}f\n" ${prec}`

# 旋转矩阵逆矩阵和平移矢量坐标的存放
for ((j=1;j<=4;j++))
do
  k=0
  for i in `awk 'NR=="'$j'"{print $0}' rotr.in`
  do
  let k++
  eval p$j['$'k]='$'i
  done
done

k=0
for temp in `awk 'NR==7{print $0}' ${file_name}`
do
  let k++
  type_num[$k]=${temp}  # 每种原子的个数
done

atom_num=`awk 'NR==7{a=0;for(i=1;i<=NF;i++){a=a+$i};print a}' ${file_name}` # 原子总数
begin=`awk 'NR>7{j=toupper(substr($1,1,1));if(j=="D"||j=="C"){print NR+1}}' ${file_name}` #开始处理的行号
let end=begin+atom_num-1  # 结束处理的行号
##################################################################
##################################################################

##################################################################
##########################   Output part   #######################
##################################################################
# 前*begin-1*行直接复制
let temp=begin-1
head -${temp} ${file_name}>>${output}
# 从*begin*行开始后面的行进行数据处理
# 首先进行坐标平移，逐行处理
# 接着对该坐标对1取余数
# 接着判断该坐标是否和1在精度prec范围内比较接近，如果是，那么将其改成0
# 然后乘以旋转矩阵的逆矩阵的对应行，然后重复上述操作

for k in $(seq ${begin} ${end})
do
  # 先平移
  temp1=`awk 'NR=="'$k'"{print $1}' ${file_name}`
  temp1=`echo "(${temp1}-(${p4[1]}))%1"|bc` # 添加平移，接着取余数
  temp=${temp1:0:1}
  if [ "${temp}" == "-" ];then
    temp1=`echo "(1+${temp1})"|bc`
  fi
  judge=`echo "(1-${temp1})/${prec}"|bc`
  if ((judge<=1));then
    temp1=0
  fi
  temp2=`awk 'NR=="'$k'"{print $2}' ${file_name}`
  temp2=`echo "(${temp2}-(${p4[2]}))%1"|bc`  # 添加平移，接着取余数
  temp=${temp2:0:1}
  if [ "${temp}" == "-" ];then
    temp2=`echo "(1+${temp2})"|bc`
  fi
  judge=`echo "(1-${temp2})/${prec}"|bc`
  if ((judge<=1));then    
    temp2=0
  fi
  temp3=`awk 'NR=="'$k'"{print $3}' ${file_name}`
  temp3=`echo "(${temp3}-(${p4[3]}))%1"|bc`  # 添加平移，接着取余数
  temp=${temp3:0:1}
  if [ "${temp}" == "-" ];then
    temp3=`echo "(1+${temp3})"|bc`
  fi
  judge=`echo "(1-${temp3})/${prec}"|bc`
  if ((judge<=1));then 
    temp3=0
  fi
  # 再旋转 
  tempp1=`echo "${p1[1]}*${temp1}+${p1[2]}*${temp2}+${p1[3]}*${temp3}"|bc`
  tempp2=`echo "${p2[1]}*${temp1}+${p2[2]}*${temp2}+${p2[3]}*${temp3}"|bc`
  tempp3=`echo "${p3[1]}*${temp1}+${p3[2]}*${temp2}+${p3[3]}*${temp3}"|bc`
# 接着调回0～1区间
  temp1=`echo "${tempp1}%1"|bc` # 取余数
  temp=${temp1:0:1}
  if [ "${temp}" == "-" ];then
    temp1=`echo "(1+${temp1})"|bc`
  fi
  judge=`echo "(1-${temp1})/${prec}"|bc`
  if ((judge<=1));then
    temp1=0
  fi
  temp2=`echo "${tempp2}%1"|bc`  # 取余数
  temp=${temp2:0:1}
  if [ "${temp}" == "-" ];then
    temp2=`echo "(1+${temp2})"|bc`
  fi
  judge=`echo "(1-${temp2})/${prec}"|bc`
  if ((judge<=1));then    
    temp2=0
  fi
  temp3=`echo "${tempp3}%1"|bc`  # 取余数
  temp=${temp3:0:1}
  if [ "${temp}" == "-" ];then
    temp3=`echo "(1+${temp3})"|bc`
  fi
  judge=`echo "(1-${temp3})/${prec}"|bc`
  if ((judge<=1));then 
    temp3=0
  fi  
  temp=`awk 'NR==8{j=toupper(substr($1,1,1));print (j=="S")?(1):(-1)}' ${file_name}` #判断是否为Selective Dynamics 模式
  if [[ "${temp}" == "1" ]];then
    temp4=`awk 'NR=="'$k'"{print $4}' ${file_name}`
    temp5=`awk 'NR=="'$k'"{print $5}' ${file_name}`
    temp6=`awk 'NR=="'$k'"{print $6}' ${file_name}`
    printf "%25.18f%25.18f%25.18f%2s%2s%2s\n" ${temp1} ${temp2} ${temp3} ${temp4} ${temp5} ${temp6}>>${output}
  else
    printf "%25.18f%25.18f%25.18f\n" ${temp1} ${temp2} ${temp3}>>${output}
  fi
done

echo "Successfully generated NEW POSCAR file!"
echo "Finished."
echo "Have a good day!"
##################################################################
##################################################################
exit 0
```



