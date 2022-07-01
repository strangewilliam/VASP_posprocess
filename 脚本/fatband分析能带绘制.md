fatband分析能带绘制

```shell
#! /bin/bash
# fat_plot
# This is a script to generate the gnuplot scipt to plot fatbands.

##################################################################
##########################  Input part   #########################
##################################################################
echo
echo "Welcome!"
echo "This is a script to generate the gnuplot scipt to plot fatbands. You need to provide:"
echo "==========================================="
echo -e "1. PBAND*.dat\n2. KLABELS"
echo "These files are generated from VASPKIT."
echo "==========================================="
echo
hom=`pwd`
if [[ "$1" != "" ]];then
  file_name=$1
else
  echo
  echo "Current work directory :" ${hom}
  echo "==========================================="
  ls
  echo "==========================================="
  echo
  read -p "Which file: " temp  # 读取PBAND文件名
  if [[ "${temp}" != "" ]];then
    file_name=${temp}
  fi
fi
if [[ ! -e "${hom}/${file_name}" ]];then
  echo "No such file!"
  echo
  echo "Quit"
  echo "Have a good day!"
  exit 1
fi
# 读取非自洽和自洽计算得到的费米能级差，用于修正费米能级
read -p "Fermi energy difference (Non-self - self): " de
if [ "${de}" == "" ];then
  de=0
fi
# 读取能带图能量范围
read -p "E_min: " y_min
read -p "E_max: " y_max
# 读取原子
read -p "Atom: " atom_type
# 高对称点标记
i=0
while [ "1" == "1" ]
do
  let i++
  let j=i+1
  K_Label[$i]=`awk -v j=$j 'NR==j{print $1}' KLABELS`
  K_Coordinate[$i]=`awk -v j=$j 'NR==j{print $2}' KLABELS`
  ASCII=`printf "%d" \'${K_Label[$i]:0:1}`
  if (( ASCII < 65)) || (( ASCII > 90 ));then
    unset -v K_Label[$i] 
    unset -v K_Coordinate[$i]
    break
  fi
done
K_path_print=""
for ((i=1;i<=${#K_Label[*]};i++))
do
  if ((i==${#K_Label[*]}));then
    K_path_print="${K_path_print} \"${K_Label[$i]}\" ${K_Coordinate[$i]}"
  else
    K_path_print="${K_path_print} \"${K_Label[$i]}\" ${K_Coordinate[$i]},"
  fi
done
# 选择要画的轨道能带
echo
echo ${file_name}
echo "==========================================="
head -2 ${file_name}
echo "==========================================="
echo
read -p "Which orbital? " temp
i=0
for orbi in ${temp}
do
  let i++ 
  orbital[$i]=${orbi}  # 轨道记号
done
# 轨道对应文件中的列号标记
tempp=`head -1 ${file_name}`
for ((i=1;i<=${#orbital[*]};i++))
do
  j=0
  for temp in ${tempp}
  do
    let j++
    if [ "${orbital[$i]}" == "${temp}" ];then
      orbital_num[$i]=$j  # 轨道对应文件中的列号
    fi
  done
  if [ "${orbital_num[$i]}" = "" ];then
    echo
    echo "No such orbital called ${orbital[$i]}!"
    echo
    echo "Quit!"
    echo "Have a good day!"
    exit 1
  fi
done
##################################################################
##########################   Output part   #######################
##################################################################
echo
echo "Generating *.gnu files ..."
for ((j=1;j<=${#orbital_num[*]};j++))
do
  title="${atom_type}_${orbital[$j]}"
  output="fat_band_${title}.gnu"
cat >${output}<<- EOF
set encoding iso_8859_1
# set terminal  postscript enhanced color font "TimesNewRoman, 11" size 5, 4
set terminal  pngcairo truecolor enhanced lw 5.0 font "TimesNewRoman, 44" size 1920, 1680
#stats file_name using (column(3))
#max = STATS_max
#set cbrange [0:max] 
set palette rgbformulae 22, 13, -33
set output '${title}.png'
set border
set title "${atom_type}\\\\_${orbital[$j]}" offset 0, -0.8 font "TimesNewRoman, 54"
set style data linespoints
unset ztics
unset key
set view 0,0
set xtics font "TimesNewRoman, 44"
set xtics offset 0, 0.3
set ytics font "TimesNewRoman, 40"
set ylabel font "TimesNewRoman, 48"
set ylabel offset 1.0, 0
set xrange [0:${K_Coordinate[${#K_Coordinate[*]}]}]
set ylabel "Energy (eV)"
set yrange [${y_min}:${y_max}]
set xtics (${K_path_print})
plot ${y_min} with filledcurves y1=${y_max} lc rgb "navy", \\
'${file_name}' u (\$1):(\$2+(${de})):(\$${orbital_num[$j]}) w lines lw 1.5 lc palette, \\
0 w l dt 2 lc rgb "red", \\
EOF
  for ((i=1;i<=${#K_Label[*]};i++))
  do
    if ((i==${#K_Label[*]}));then
      echo "'< echo \"${K_Coordinate[$i]} ${y_min} \\n ${K_Coordinate[$i]} ${y_max}\"' w l dt 2 lc rgb \"gray\"">> ${output}
    else
      echo "'< echo \"${K_Coordinate[$i]} ${y_min} \\n ${K_Coordinate[$i]} ${y_max}\"' w l dt 2 lc rgb \"gray\", \\">> ${output}
    fi
  done
done
echo
echo "*.gnu files:"
echo "==========================================="
ls *.gnu
echo "==========================================="
echo
echo "Finished."
echo "Have a good day!"
echo
##################################################################
##################################################################
exit 0
```
