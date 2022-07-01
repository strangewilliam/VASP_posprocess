# POSCAR等价性判断脚本

 posim

```shell
#! /bin/bash
# posim
# This is a script to judge the smae POSCAR files provided in the selected folder.
# You need the script *posarrange+PERIODIC_TABLE*  *posorder*  *postrans* and be under Phonopy environment.
# The PERIODIC_TABLE has to be put in the *${program}* directory

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
  judge_1=`awk -v abs=${abs} 'BEGIN{print (abs<1)?0:1}'`
  abs=`awk -v p=${pre} -v abs=${abs} 'BEGIN{a=p-abs;print a}'`
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
program="/fs12/home/zhj_xujq/bin"
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
prec=1E-7
if [[ "$2" != "" ]];then
  prec=$2
else
  read -p "Precision[Default=1E-7]: " temp
  if [ "${temp}" != "" ];then
    prec=${temp}
  fi
fi
echo
rm -rf ${hom}/com
rm -rf ${hom}/workout
touch ${work}/log
touch ${work}/order.out
touch ${work}/equi.out
echo |tee -a ${work}/log
echo "Precision: ${prec}"|tee -a ${work}/log
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
# ======加速计算：后面和前面等价的POSCAR文件后续就不需要参与比较了 ====== #
no=0
NONEED=""
# ================================================================= #
for ((i=1;i<num;i++))
do
  # ======加速计算：后面和前面等价的POSCAR文件后续就不需要参与比较了 ====== #
  for skip in ${NONEED[*]}
  do
    if ((i==skip));then
      echo "No need to compare POSCAR_${i} with others"|tee -a ${work}/log
      continue 2
    fi
  done
  # ================================================================= #
  for ((j=i+1;j<=num;j++))
  do
    # ======加速计算：后面和前面等价的POSCAR文件后续就不需要参与比较了 ====== #
    for skip in ${NONEED[*]}
    do
      if ((j==skip));then
        echo "No need to compare POSCAR_${j} with others"|tee -a ${work}/log
        continue 2
      fi
    done
    # ================================================================= #
    let order=(i-1)*num+j
    judge[${order}]=0
    dir="${hom}/compare_${i}_${j}"
    mkdir ${dir}
    cp -pai ${hom}/${file[$i]} ${dir}/POSCAR_${i}
    cp -pai ${hom}/${file[$j]} ${dir}/POSCAR_${j}
    file1_p=`printf "%16s" ${file[$i]}`
    file2_p=`printf "%16s" ${file[$j]}`
    echo|tee -a ${work}/log
    echo "compare_${i}_${j}"|tee -a ${work}/log
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
      rm -rf POSCAR_g_${k}_i
      posarrange POSCAR_g_${k} POSCAR_g_${k}_i p ${program} >>${work}/log
      mv POSCAR_g_${k}_i POSCAR_g_${k}
      phonopy -c POSCAR_g_${k} --symmetry>sym_${k}   #注意：此时的sym_*文件是基于POSCAR_g_*的
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
        same ${num1} ${num2} ${prec}
        temp=$?
        let judge[${order}]=judge[${order}]+temp
        if (( judge[${order}] != 0 ));then
          echo "Phonopy转化生成的二者原胞的基矢不同!"        
          continue 3
        fi
      done
    done
    echo "Phonopy转化生成的原胞的基矢一致!" 
    # 比较一下原子总数是否一致 
    temp=`awk 'ARGIND==1{if(FNR==7){a=0;for(i=1;i<=NF;i++){a=a+$i}}}ARGIND==2{if(FNR==7){b=0;for(i=1;i<=NF;i++){b=b+$i}}}END{print (a!=b)?1:0}' ${file1} ${file2}`
    let judge[${order}]=judge[${order}]+temp
    if (( judge[${order}] != 0 ));then
      continue
    fi
    echo "原子总数一致!"    
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
          eval printf "%3s%4s%3s%3d'\n'" '$'{atom${type}_order['$'k]} '$'{atom_type['$'{type}]} '$'{atom${type}_Wyckoff_mark['$'k]} '$'{atom${type}_equiv_num['$'k]}>> ${dir}/sym_log_${ii}
          eval printf "%3s%4s%3s%3d%3d'\n'" '$'{atom${type}_order['$'k]} '$'{atom_type['$'{type}]} '$'{atom${type}_Wyckoff_mark['$'k]} '$'{atom${type}_equiv_num['$'k]} '$'{atom${type}_irre['$'k]}>> ${dir}/sym_log_d_${ii}
        done
      done
    done
    echo "对称性分析文件sym*已经生成！"
    # 比较生成的两个sym_log文件是不是一样的
    file1=${dir}/sym_log_${i}
    file2=${dir}/sym_log_${j}
    # 先比较行数
    lines1=`cat ${file1} | wc -l`
    lines2=`cat ${file2} | wc -l`
    if (( lines1 != lines2 ));then
      let judge[${order}]++
      echo "二者的sym_log文件行数不一样!"
      continue
    fi
    echo "二者的sym_log文件行数一样!"
    # 再比较细节
    for ((k=1;k<=lines1;k++))
    do  
      line1=`awk 'NR=="'$k'"{print $0}' ${file1}`
      line2=`awk 'NR=="'$k'"{print $0}' ${file2}` 
      if [[ "${line1}" != "${line2}" ]];then
        let judge[${order}]++
        echo "二者的sym_log文件细节不一样!"
        continue 2 
      fi    
    done 
    echo "二者的sym_log文件是一样的!"
    # 以上比较了PPOSCAR的晶格常数、原子总数、对称性信息（包含了原子种类及其个数）。如果说都是一样的话，那么接下来就要具体比较POSCAR文件了。策略：把POSCAR文件中首个元素中字母最小、重数最低的那类原子（对应于sym_log文件中的第一个原子）放到原点位置，比较两个POSCAR文件。如果遇到不等价原子数目一致的情况，按照原先POSCAR文件中的元素排列顺序来（要求两个POSCAR文件中元素排列顺序一致），如果遇到字母、重数都一样的情况，那么就固定其中的一个POSCAR，让另一个POSCAR文件把这些看起来信息一模一样的不等价位置都位移至原点尝试一下。如果还是不一样，就变动第一个POSCAR文件的原子种类，重复以上步骤。
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
    begin=`awk 'NR>7{j=toupper(substr($1,1,1));if(j=="D"||j=="C"){print NR+1}}' POSCAR_g_${i}`  #开始处理的行号    
    let end=begin+total_num-1  # 结束处理的行号
    
    num_select_temp=0
    begin_add=${begin}
    for type_i in $(seq 1 ${#atom_type[*]})
    do 
      atom_origin=`awk -v type_i=${atom_type[${type_i}]} '{if($2==type_i){print $5;exit}}' ${dir}/sym_log_d_$i` # 确定要放到原点的原子对应的序号
      let temp=begin+atom_origin-1 #要放到原点的原子对应的行号
      # 读取原点原子位置信息
      x_move=`awk 'NR=="'${temp}'"{print $1}' POSCAR_g_${i}`
      y_move=`awk 'NR=="'${temp}'"{print $2}' POSCAR_g_${i}`
      z_move=`awk 'NR=="'${temp}'"{print $3}' POSCAR_g_${i}`
#echo begin+atom_origin-1
#echo $temp
#echo i= $i
#echo atom_type: ${atom_type[${type_i}]}
#echo x_move: $x_move
#echo y_move: $y_move
#echo z_move: $z_move
      # 把POSCAR_g_i文件移动选定原子到原点
      rm -rf POSCAR_g_${i}_t*
      postrans POSCAR_g_${i} POSCAR_g_${i}_t ${prec} ${x_move} ${y_move} ${z_move}>>${work}/log
      # POSCAR 文件排序
      posorder POSCAR_g_${i}_t POSCAR_g_${i}_t_o ${prec}>>${work}/log
      judge_temp=0
      let begin_add=begin_add+num_select_temp
      num_select_temp=`awk -v i=${type_i} 'NR==7{print $i}' POSCAR_g_${i}`
      for ((tempp=begin_add;tempp<=begin_add+num_select_temp-1;tempp++))
      do
        # 读取原点原子位置信息
        x_move=`awk 'NR=="'${tempp}'"{print $1}' POSCAR_g_${j}`
        y_move=`awk 'NR=="'${tempp}'"{print $2}' POSCAR_g_${j}`
        z_move=`awk 'NR=="'${tempp}'"{print $3}' POSCAR_g_${j}`
#echo begin_add+num_select-1:
#let a=begin_add+num_select-1
#echo $a
#echo num_select: $num_select
#echo tempp: $tempp
#echo j= $j
#echo x_move: $x_move
#echo y_move: $y_move
#echo z_move: $z_move       
        # 把POSCAR_g_j文件移动选定原子到原点
        rm -rf POSCAR_g_${j}_t*
        postrans POSCAR_g_${j} POSCAR_g_${j}_t ${prec} ${x_move} ${y_move} ${z_move}>>${work}/log
        # POSCAR 文件排序
        posorder POSCAR_g_${j}_t POSCAR_g_${j}_t_o ${prec}>>${work}/log
        # 比较一下POSCAR文件每个原子位置是不是完全一样的
        file1="POSCAR_g_${i}_t_o"
        file2="POSCAR_g_${j}_t_o"
        num_true=`awk 'END{print NR}' ${file1}`
        for ((k1=9;k1<=${num_true};k1++))
        do
          for l in 1 2 3
          do
            num1=`awk -v l=${l} 'NR=="'${k1}'"{print $l}' ${file1}`
            num2=`awk -v l=${l} 'NR=="'${k1}'"{print $l}' ${file2}`
            same ${num1} ${num2} ${prec}
            temp=$?
            let judge_temp=judge_temp+temp
            if [ "${judge_temp}" != "0" ];then        
              break 2
            fi
          done
        done
        if [ "${judge_temp}" != "0" ];then
          judge_temp=0
          continue
        else
          judge_temp=1
          break 2
        fi
      done  
    done
    if [ "${judge_temp}" == "1" ];then
      judge[${order}]=eq
      echo "POSCAR文件等价！"
      # ======加速计算：后面和前面等价的POSCAR文件后续就不需要参与比较了 ====== #
      let no++
      NONEED[$no]=$j
      # ================================================================= #
    else
      judge[${order}]=1
      echo "POSCAR文件不等价！"
    fi
    cd ${hom}   
  done
done
##################################################################
##################################################################

##################################################################
##########################   Output part   #######################
##################################################################
mkdir ${hom}/com
echo
echo "Moving *compare_** to *com* ..."
mv ${hom}/compare_* ${hom}/com
echo
echo|tee -a ${work}/log
echo "Equivalent file"|tee -a ${work}/log
echo "==========================================="|tee -a ${work}/log
for ((i=1;i<num;i++))
do
  for ((j=i+1;j<=num;j++))
  do
    let order=(i-1)*num+j
    if [ "${judge[${order}]}" == "eq" ];then
      printf "%16s%16s\n" ${file[$i]} ${file[$j]}| tee -a ${work}/equi.out
    fi
  done
done
cat ${work}/equi.out >>${work}/log
echo "==========================================="|tee -a ${work}/log
echo|tee -a ${work}/log
mv ${work} ${hom}
work=${hom}/workout
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





posimp

```shell
#! /bin/bash
# posimp
# This is a script to judge the smae POSCAR files provided in the selected folder.
# The POSCAR files may be different in the positions of some type of atoms and the number of each type of atom are the same.
# You need the script *posarrange+PERIODIC_TABLE*  *posorder*  *postrans* and be under Phonopy environment.
# The PERIODIC_TABLE has to be put in the *${program}* directory

##################################################################
########################  function part   ########################
##################################################################
function same(){
  local ppre=$3
  local ppre=${ppre:0-1:1}
  local abs=`echo "scale=${ppre};($1-$2)/1"|bc`
  local temp=${abs:0:1}
  if [ ${temp} == "-" ];then
    abs=`echo "0-(${abs})"|bc`
  fi
  local judge=`awk -v abs=${abs} 'BEGIN{print (abs==0)?(0):(1)}'`
  return ${judge}
}
##################################################################
##################################################################

##################################################################
##########################  Input part   #########################
##################################################################
hom=`pwd`
program="/fs12/home/zhj_xujq/bin"
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
prec=1E-7
if [[ "$2" != "" ]];then
  prec=$2
else
  read -p "Precision[Default=1E-7]: " temp
  if [ "${temp}" != "" ];then
    prec=${temp}
  fi
fi

# 输入需要比较的元素
if [[ "$3" != "" ]];then
  i=0
  while [ "1" == "1" ]
  do
    let i++
    let i_temp=2+i
    eval  temp='$'${i_temp}
    if [ "${temp}" != "" ];then
      atom_type_selected[$i]=${temp} # 需要比较的元素符号
    else
      break
    fi
  done
else
  read -p "The elements needed to be compare: " temp
  if [ "${temp}" != "" ];then
    i=0
    for things in ${temp}
    do
      let i++
      atom_type_selected[$i]=${things} # 需要比较的元素符号 
    done
  else
    echo "No atom input!"
    echo
    echo "Quit"
    echo "Have a good day!"
    exit 1
  fi
fi

echo
rm -rf ${hom}/com*
rm -rf ${hom}/workout
touch ${work}/log
touch ${work}/order.out
touch ${work}/equi.out
echo |tee -a ${work}/log
echo "Precision: ${prec}"|tee -a ${work}/log
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
# ======加速计算：后面和前面等价的POSCAR文件后续就不需要参与比较了 ====== #
no=0
NONEED=""
# ================================================================= #
for ((i=1;i<num;i++))
do
  # ======加速计算：后面和前面等价的POSCAR文件后续就不需要参与比较了 ====== #
  for skip in ${NONEED[*]}
  do
    if ((i==skip));then
      echo "No need to compare POSCAR_${i} with others"|tee -a ${work}/log
      continue 2
    fi
  done
  # ================================================================= #
  for ((j=i+1;j<=num;j++))
  do
    # ======加速计算：后面和前面等价的POSCAR文件后续就不需要参与比较了 ====== #
    for skip in ${NONEED[*]}
    do
      if ((j==skip));then
        echo "No need to compare POSCAR_${j} with others"|tee -a ${work}/log
        continue 2
      fi
    done
    # ================================================================= #
    let order=(i-1)*num+j
    judge[${order}]=0
    dir="${hom}/compare_${i}_${j}"
    mkdir ${dir}
    cp -pai ${hom}/${file[$i]} ${dir}/POSCAR_${i}
    cp -pai ${hom}/${file[$j]} ${dir}/POSCAR_${j}
    file1_p=`printf "%16s" ${file[$i]}`
    file2_p=`printf "%16s" ${file[$j]}`
    echo|tee -a ${work}/log
    echo "compare_${i}_${j}"|tee -a ${work}/log
    echo "==========================================="|tee -a ${work}/log
    echo "${file1_p} =======>  POSCAR_${i}"|tee -a ${work}/log
    echo "${file2_p} =======>  POSCAR_${j}"|tee -a ${work}/log
    echo "==========================================="|tee -a ${work}/log
    echo|tee -a ${work}/log
    cd ${dir}
    file1="POSCAR_${i}"
    file2="POSCAR_${j}"    
    posinit ${file1} ${file1}_i >>${work}/log
    posinit ${file2} ${file2}_i >>${work}/log
    posarrange ${file1}_i ${file1}_r p ${program} >>${work}/log
    posarrange ${file2}_i ${file2}_r p ${program} >>${work}/log
    # 接下来就使用 *_r 文件了
    file1="${file1}_r"
    file2="${file2}_r"
    # 比较第6、7行
    for ((k=6;k<=7;k++))
    do  
      line1=`awk 'NR=="'$k'"{print $0}' ${file1}`
      line2=`awk 'NR=="'$k'"{print $0}' ${file2}` 
      if [[ "${line1}" != "${line2}" ]];then
        let judge[${order}]++
        echo "二者的元素种类或者原子数目不一样!"|tee -a ${work}/log
        continue 2 
      fi    
    done 
    echo "二者的元素种类和原子数目一致!"|tee -a ${work}/log
    
    k=0;i1=0
    for temp in `awk 'NR==6{print $0}' ${file1}`
    do
      let k++
      atom_type[$k]=${temp}  # 原子种类
      for type in ${atom_type_selected[*]}
      do
        if [ "${type}" == "${temp}" ];then
          let i1++
          atom_selected[${i1}]=$k # 需要比较的元素在 *_r 中的位置序号
          break
        fi
      done
    done
    # 判断选定原子格式是否正确
    if ((${#atom_selected[*]}!=${#atom_type_selected[*]}));then
      echo "Wrong input format of the atom selected!"
      echo
      echo "Quit"
      echo "Have a good day!"
      exit 1
    fi
    k=0
    total_num=0 # 原子总数
    for temp in `awk 'NR==7{print $0}' ${file1}`
    do
      let k++
      type_num[$k]=${temp}  # 原子数目
      let total_num=total_num+temp
    done
    # 判断该方法是否适用
    if ((total_num<4));then
      echo "WARNNING: the total number of atoms is less than 4. The results may be incorrect!"|tee -a ${work}/log
      echo
    fi    
    temp=1
    for ((i1=1;i1<=${#type_num[*]};i1++))
    do
      atom_begin[$i1]=${temp} # 每个元素起始序号
      let atom_end[$i1]=atom_begin[$i1]+type_num[$i1]-1  # 每个元素终止序号
      let temp=atom_end[$i1]+1
    done
    # 分元素种类进行比较
    for ((i1=1;i1<=${#atom_selected[*]};i1++))
    do
      for ((j1=atom_begin[${atom_selected[$i1]}];j1<=atom_end[${atom_selected[$i1]}];j1++)) # j1 --需要比较的原子序号
      do
        name="${atom_type[${atom_selected[$i1]}]}_${j1}"
        posdist ${file1} ${j1} ${file1}_dist_${name} ${prec}>>${work}/log
        posdist ${file2} ${j1} ${file2}_dist_${name} ${prec}>>${work}/log 
        file1_dist="${file1}_dist_${name}"
        file2_dist="${file2}_dist_${name}"
        # 把距离从文件中读出来
        for ((k=1;k<=total_num;k++))
        do
          let k_temp=k+1
          temp=`awk -v k=${k_temp} 'NR==k{print $2}' ${file1_dist}`
          eval file1_${name}['$'k]='$'{temp}
          temp=`awk -v k=${k_temp} 'NR==k{print $2}' ${file2_dist}`
          eval file2_${name}['$'k]='$'{temp}
        done
        # 对这两个数组排序
        for ((k0=1;k0<=${#type_num[*]};k0++))
        do
          for ((k1=atom_begin[${k0}];k1<=atom_end[${k0}]-1;k1++))
          do
            for ((k2=k1+1;k2<=atom_end[${k0}];k2++))
            do
              eval w1='$'{file1_${name}['$'k1]}
              eval w2='$'{file1_${name}['$'k2]} 
              judge=`awk -v w1=${w1} -v w2=${w2} 'BEGIN{print (w1>w2)?(1):(-1)}'`
              if ((judge>0));then
                eval temp='$'{file1_${name}['$'k1]}
                eval file1_${name}['$'k1]='$'{file1_${name}['$'k2]}
                eval file1_${name}['$'k2]='$'{temp}
              fi
              eval w1='$'{file2_${name}['$'k1]}
              eval w2='$'{file2_${name}['$'k2]}              
              judge=`awk -v w1=${w1} -v w2=${w2} 'BEGIN{print (w1>w2)?(1):(-1)}'`
              if ((judge>0));then
                eval temp='$'{file2_${name}['$'k1]}
                eval file2_${name}['$'k1]='$'{file2_${name}['$'k2]}
                eval file2_${name}['$'k2]='$'{temp}
              fi
            done
          done
        done      
      done
      # 比较两个文件中的数组  
      name="${atom_type[${atom_selected[$i1]}]}"
      unset NONEED_1;no_1=0
      for ((j1=atom_begin[${atom_selected[$i1]}];j1<=atom_end[${atom_selected[$i1]}];j1++))
      do
        judge_temp=0
        for ((j2=atom_begin[${atom_selected[$i1]}];j2<=atom_end[${atom_selected[$i1]}];j2++))
        do 
         # ========================加速计算========================== #
          for skip_1 in ${NONEED_1[*]}
          do
            if ((j2==skip_1));then
              continue 2
            fi
          done
          # ========================================================== #
          for ((k1=1;k1<=total_num;k1++))
          do
            eval num1='$'{file1_${name}_${j1}['$'k1]}
            eval num2='$'{file2_${name}_${j2}['$'k1]}
           same ${num1} ${num2} ${prec}
            temp=$?
            if [ "${temp}" != "0" ];then   
              let judge_temp++
              break
            fi
          done
          if [ "${judge_temp}" == "0" ];then 
            # ========================加速计算========================== #
            let no_1++
            NONEED_1[${no_1}]=$j2
            # ========================================================== #
            break
          fi
        done
        if ((judge_temp>=type_num[${atom_selected[$i1]}]));then
          judge_temp=0
          break 2
        else
          judge_temp=1
        fi
      done     
    done  
    if [ "${judge_temp}" == "1" ];then
      judge[${order}]=eq
      echo "POSCAR文件等价!"|tee -a ${work}/log
      # ======加速计算：后面和前面等价的POSCAR文件后续就不需要参与比较了 ====== #
      let no++
      NONEED[$no]=$j
      # ================================================================= #
    else
      judge[${order}]=1
      echo "POSCAR文件不等价！"
    fi
    cd ${hom}   
  done
done
##################################################################
##################################################################

##################################################################
##########################   Output part   #######################
##################################################################
mkdir ${hom}/com
echo
echo "Moving *compare_** to *com* ..."
mv ${hom}/compare_* ${hom}/com
echo
echo|tee -a ${work}/log
echo "Equivalent file"|tee -a ${work}/log
echo "==========================================="|tee -a ${work}/log
for ((i=1;i<num;i++))
do
  for ((j=i+1;j<=num;j++))
  do
    let order=(i-1)*num+j
    if [ "${judge[${order}]}" == "eq" ];then
      printf "%16s%16s\n" ${file[$i]} ${file[$j]}| tee -a ${work}/equi.out
    fi
  done
done
cat ${work}/equi.out >>${work}/log
echo "==========================================="|tee -a ${work}/log
echo|tee -a ${work}/log
mv ${work} ${hom}
work=${hom}/workout
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



```shell
#! /bin/bash
# poscom
# This is a script to judge the smae POSCAR files provided in the selected folder.
# The POSCAR files may be different in the positions of some type of atoms and the number of each type of atom are the same.
# You need the script *posarrange+PERIODIC_TABLE*  *posorder*  *postrans* and be under Phonopy environment.
# The PERIODIC_TABLE has to be put in the *${program}* directory

##################################################################
########################  function part   ########################
##################################################################
function same(){
  local ppre=$3
  local ppre=${ppre:0-1:1}
  local abs=`echo "scale=${ppre};($1-$2)/1"|bc`
  local temp=${abs:0:1}
  if [ ${temp} == "-" ];then
    abs=`echo "0-(${abs})"|bc`
  fi
  local judge=`awk -v abs=${abs} 'BEGIN{print (abs==0)?(0):(1)}'`
  return ${judge}
}

function same1(){
  local ppre=$3
  local ppre=${ppre:0-1:1}
  local pre=1
  for ((fj=1;fj<=${ppre};fj++))
  do
    pre="${pre}0"
  done   
  local abs=`awk -v p=${pre} -v n1=$1 -v n2=$2 'BEGIN{a=p*(n1-n2);print a}'` 
  abs=`awk -v abs=${abs} 'BEGIN{print (abs<0)?(-abs):(abs)}'`
  local judge_1=`awk -v abs=${abs} 'BEGIN{print (abs<1)?0:1}'`
  abs=`awk -v p=${pre} -v abs=${abs} 'BEGIN{a=p-abs;print a}'`
  abs=`awk -v abs=${abs} 'BEGIN{print (abs<0)?(-abs):(abs)}'`  
  local judge_2=`awk -v abs=${abs} 'BEGIN{print (abs<1)?0:1}'`
  let judge=(judge_1+judge_2)/2
  return ${judge}
}
##################################################################
##################################################################

##################################################################
##########################  Input part   #########################
##################################################################
hom=`pwd`
program="/fs12/home/zhj_xujq/bin"
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
prec=1E-7
if [[ "$2" != "" ]];then
  prec=$2
else
  read -p "Precision[Default=1E-7]: " temp
  if [ "${temp}" != "" ];then
    prec=${temp}
  fi
fi

# 输入需要比较的元素
if [[ "$3" != "" ]];then
  i=0
  while [ "1" == "1" ]
  do
    let i++
    let i_temp=2+i
    eval  temp='$'${i_temp}
    if [ "${temp}" != "" ];then
      atom_type_selected[$i]=${temp} # 需要比较的元素符号
    else
      break
    fi
  done
else
  read -p "The elements needed to be compare: " temp
  if [ "${temp}" != "" ];then
    i=0
    for things in ${temp}
    do
      let i++
      atom_type_selected[$i]=${things} # 需要比较的元素符号 
    done
  else
    echo "No atom input!"
    echo
    echo "Quit"
    echo "Have a good day!"
    exit 1
  fi
fi

echo
rm -rf ${hom}/com*
rm -rf ${hom}/workout
touch ${work}/log
touch ${work}/order.out
touch ${work}/equi.out
echo |tee -a ${work}/log
echo "Precision: ${prec}"|tee -a ${work}/log
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
# ======加速计算：后面和前面等价的POSCAR文件后续就不需要参与比较了 ====== #
no=0
NONEED=""
# ================================================================= #
for ((i=1;i<num;i++))
do
  # ======加速计算：后面和前面等价的POSCAR文件后续就不需要参与比较了 ====== #
  for skip in ${NONEED[*]}
  do
    if ((i==skip));then
      echo "No need to compare POSCAR_${i} with others"|tee -a ${work}/log
      continue 2
    fi
  done
  # ================================================================= #
  for ((j=i+1;j<=num;j++))
  do
    # ======加速计算：后面和前面等价的POSCAR文件后续就不需要参与比较了 ====== #
    for skip in ${NONEED[*]}
    do
      if ((j==skip));then
        echo "No need to compare POSCAR_${j} with others"|tee -a ${work}/log
        continue 2
      fi
    done
    # ================================================================= #
    let order=(i-1)*num+j
    judge[${order}]=0
    dir="${hom}/compare_${i}_${j}"
    mkdir ${dir}
    cp -pai ${hom}/${file[$i]} ${dir}/POSCAR_${i}
    cp -pai ${hom}/${file[$j]} ${dir}/POSCAR_${j}
    file1_p=`printf "%16s" ${file[$i]}`
    file2_p=`printf "%16s" ${file[$j]}`
    echo|tee -a ${work}/log
    echo "compare_${i}_${j}"|tee -a ${work}/log
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
      rm -rf POSCAR_g_${k}_i
      posarrange POSCAR_g_${k} POSCAR_g_${k}_i p ${program} >>${work}/log
      mv POSCAR_g_${k}_i POSCAR_g_${k}
      phonopy -c POSCAR_g_${k} --symmetry>sym_${k}   #注意：此时的sym_*文件是基于POSCAR_g_*的
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
        same1 ${num1} ${num2} ${prec}
        temp=$?
        let judge[${order}]=judge[${order}]+temp
        if (( judge[${order}] != 0 ));then
          echo "Phonopy转化生成的二者原胞的基矢不同!"        
          continue 3
        fi
      done
    done
    echo "Phonopy转化生成的原胞的基矢一致!" 
    # 比较一下原子总数是否一致 
    temp=`awk 'ARGIND==1{if(FNR==7){a=0;for(i=1;i<=NF;i++){a=a+$i}}}ARGIND==2{if(FNR==7){b=0;for(i=1;i<=NF;i++){b=b+$i}}}END{print (a!=b)?1:0}' ${file1} ${file2}`
    let judge[${order}]=judge[${order}]+temp
    if (( judge[${order}] != 0 ));then
      continue
    fi
    echo "原子总数一致!"    
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
        type_num[$k]=${temp}  # 原子数目
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
      for ((type=1;type<=${#type_num[*]};type++))
      do
        let limit_up=limit_bot+type_num[${type}]-1
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
      for ((type=1;type<=${#type_num[*]};type++))
      do
        eval declare -a atom${type}_order
        for ((k=1;k<=type_num[${type}];k++))
        do 
          eval atom${type}_order['$'k]='$'k
        done
      done
      # 排序：先按照Wyckoff记号，a b c d ..., 然后按照每个Wyckoff位置含有的原子的数目从低到高排列
      for ((type=1;type<=${#type_num[*]};type++))
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
      for ((type=1;type<=${#type_num[*]};type++))
      do
        for ((k=1;k<=type_num[${type}];k++))
        do      
          eval printf "%3s%4s%3s%3d'\n'" '$'{atom${type}_order['$'k]} '$'{atom_type['$'{type}]} '$'{atom${type}_Wyckoff_mark['$'k]} '$'{atom${type}_equiv_num['$'k]}>> ${dir}/sym_log_${ii}
          eval printf "%3s%4s%3s%3d%3d'\n'" '$'{atom${type}_order['$'k]} '$'{atom_type['$'{type}]} '$'{atom${type}_Wyckoff_mark['$'k]} '$'{atom${type}_equiv_num['$'k]} '$'{atom${type}_irre['$'k]}>> ${dir}/sym_log_d_${ii}
        done
      done
    done
    echo "对称性分析文件sym*已经生成！"
    # 比较生成的两个sym_log文件是不是一样的
    file1=${dir}/sym_log_${i}
    file2=${dir}/sym_log_${j}
    # 先比较行数
    lines1=`cat ${file1} | wc -l`
    lines2=`cat ${file2} | wc -l`
    if (( lines1 != lines2 ));then
      let judge[${order}]++
      echo "二者的sym_log文件行数不一样!"
      continue
    fi
    echo "二者的sym_log文件行数一样!"
    # 再比较细节
    for ((k=1;k<=lines1;k++))
    do  
      line1=`awk 'NR=="'$k'"{print $0}' ${file1}`
      line2=`awk 'NR=="'$k'"{print $0}' ${file2}` 
      if [[ "${line1}" != "${line2}" ]];then
        let judge[${order}]++
        echo "二者的sym_log文件细节不一样!"
        continue 2 
      fi    
    done 
    echo "二者的sym_log文件是一样的!"
    
    file1="POSCAR_${i}"
    file2="POSCAR_${j}"    
    posinit ${file1} ${file1}_i >>${work}/log
    posinit ${file2} ${file2}_i >>${work}/log
    posarrange ${file1}_i ${file1}_r p ${program} >>${work}/log
    posarrange ${file2}_i ${file2}_r p ${program} >>${work}/log
    # 接下来就使用 *_r 文件了
    file1="${file1}_r"
    file2="${file2}_r"
    # 比较第6、7行
    for ((k=6;k<=7;k++))
    do  
      line1=`awk 'NR=="'$k'"{print $0}' ${file1}`
      line2=`awk 'NR=="'$k'"{print $0}' ${file2}` 
      if [[ "${line1}" != "${line2}" ]];then
        let judge[${order}]++
        echo "二者的元素种类或者原子数目不一样!"|tee -a ${work}/log
        continue 2 
      fi    
    done 
    echo "二者的元素种类和原子数目一致!"|tee -a ${work}/log
    
    k=0;i1=0
    for temp in `awk 'NR==6{print $0}' ${file1}`
    do
      let k++
      atom_type[$k]=${temp}  # 原子种类
      for type in ${atom_type_selected[*]}
      do
        if [ "${type}" == "${temp}" ];then
          let i1++
          atom_selected[${i1}]=$k # 需要比较的元素在 *_r 中的位置序号
          break
        fi
      done
    done
    # 判断选定原子格式是否正确
    if ((${#atom_selected[*]}!=${#atom_type_selected[*]}));then
      echo "Wrong input format of the atom selected!"
      echo
      echo "Quit"
      echo "Have a good day!"
      exit 1
    fi
    k=0
    total_num=0 # 原子总数
    for temp in `awk 'NR==7{print $0}' ${file1}`
    do
      let k++
      type_num[$k]=${temp}  # 原子数目
      let total_num=total_num+temp
    done
    # 判断该方法是否适用
    if ((total_num<4));then
      echo "WARNNING: the total number of atoms is less than 4. The results may be incorrect!"|tee -a ${work}/log
      echo
    fi    
    temp=1
    for ((i1=1;i1<=${#type_num[*]};i1++))
    do
      atom_begin[$i1]=${temp} # 每个元素起始序号
      let atom_end[$i1]=atom_begin[$i1]+type_num[$i1]-1  # 每个元素终止序号
      let temp=atom_end[$i1]+1
    done
    # 分元素种类进行比较
    for ((i1=1;i1<=${#atom_selected[*]};i1++))
    do
      for ((j1=atom_begin[${atom_selected[$i1]}];j1<=atom_end[${atom_selected[$i1]}];j1++)) # j1 --需要比较的原子序号
      do
        name="${atom_type[${atom_selected[$i1]}]}_${j1}"
        posdist ${file1} ${j1} ${file1}_dist_${name} ${prec}>>${work}/log
        posdist ${file2} ${j1} ${file2}_dist_${name} ${prec}>>${work}/log 
        file1_dist="${file1}_dist_${name}"
        file2_dist="${file2}_dist_${name}"
        # 把距离从文件中读出来
        for ((k=1;k<=total_num;k++))
        do
          let k_temp=k+1
          temp=`awk -v k=${k_temp} 'NR==k{print $2}' ${file1_dist}`
          eval file1_${name}['$'k]='$'{temp}
          temp=`awk -v k=${k_temp} 'NR==k{print $2}' ${file2_dist}`
          eval file2_${name}['$'k]='$'{temp}
        done
        # 对这两个数组排序
        for ((k0=1;k0<=${#type_num[*]};k0++))
        do
          for ((k1=atom_begin[${k0}];k1<=atom_end[${k0}]-1;k1++))
          do
            for ((k2=k1+1;k2<=atom_end[${k0}];k2++))
            do
              eval w1='$'{file1_${name}['$'k1]}
              eval w2='$'{file1_${name}['$'k2]} 
              judge=`awk -v w1=${w1} -v w2=${w2} 'BEGIN{print (w1>w2)?(1):(-1)}'`
              if ((judge>0));then
                eval temp='$'{file1_${name}['$'k1]}
                eval file1_${name}['$'k1]='$'{file1_${name}['$'k2]}
                eval file1_${name}['$'k2]='$'{temp}
              fi
              eval w1='$'{file2_${name}['$'k1]}
              eval w2='$'{file2_${name}['$'k2]}              
              judge=`awk -v w1=${w1} -v w2=${w2} 'BEGIN{print (w1>w2)?(1):(-1)}'`
              if ((judge>0));then
                eval temp='$'{file2_${name}['$'k1]}
                eval file2_${name}['$'k1]='$'{file2_${name}['$'k2]}
                eval file2_${name}['$'k2]='$'{temp}
              fi
            done
          done
        done      
      done
      # 比较两个文件中的数组  
      name="${atom_type[${atom_selected[$i1]}]}"
      unset NONEED_1;no_1=0
      for ((j1=atom_begin[${atom_selected[$i1]}];j1<=atom_end[${atom_selected[$i1]}];j1++))
      do
        judge_temp=0
        for ((j2=atom_begin[${atom_selected[$i1]}];j2<=atom_end[${atom_selected[$i1]}];j2++))
        do 
         # ========================加速计算========================== #
          for skip_1 in ${NONEED_1[*]}
          do
            if ((j2==skip_1));then
              continue 2
            fi
          done
          # ========================================================== #
          for ((k1=1;k1<=total_num;k1++))
          do
            eval num1='$'{file1_${name}_${j1}['$'k1]}
            eval num2='$'{file2_${name}_${j2}['$'k1]}
           same ${num1} ${num2} ${prec}
            temp=$?
            if [ "${temp}" != "0" ];then   
              let judge_temp++
              break
            fi
          done
          if [ "${judge_temp}" == "0" ];then 
            # ========================加速计算========================== #
            let no_1++
            NONEED_1[${no_1}]=$j2
            # ========================================================== #
            break
          fi
        done
        if ((judge_temp>=type_num[${atom_selected[$i1]}]));then
          judge_temp=0
          break 2
        else
          judge_temp=1
        fi
      done     
    done  
    if [ "${judge_temp}" == "1" ];then
      judge[${order}]=eq
      echo "POSCAR文件等价!"|tee -a ${work}/log
      # ======加速计算：后面和前面等价的POSCAR文件后续就不需要参与比较了 ====== #
      let no++
      NONEED[$no]=$j
      # ================================================================= #
    else
      judge[${order}]=1
      echo "POSCAR文件不等价！"
    fi
    cd ${hom}   
  done
done
##################################################################
##################################################################

##################################################################
##########################   Output part   #######################
##################################################################
mkdir ${hom}/com
echo
echo "Moving *compare_** to *com* ..."
mv ${hom}/compare_* ${hom}/com
echo
echo|tee -a ${work}/log
echo "Equivalent file"|tee -a ${work}/log
echo "==========================================="|tee -a ${work}/log
for ((i=1;i<num;i++))
do
  for ((j=i+1;j<=num;j++))
  do
    let order=(i-1)*num+j
    if [ "${judge[${order}]}" == "eq" ];then
      printf "%16s%16s\n" ${file[$i]} ${file[$j]}| tee -a ${work}/equi.out
    fi
  done
done
cat ${work}/equi.out >>${work}/log
echo "==========================================="|tee -a ${work}/log
echo|tee -a ${work}/log
mv ${work} ${hom}
work=${hom}/workout
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

