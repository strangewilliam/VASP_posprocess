# DOS 绘制

## TDOS

```shell
#! /bin/bash
# tdos_plot
# This is a script to generate the gnuplot scipt to plot TOTAL DOS.

##################################################################
##########################  Input part   #########################
##################################################################
echo
echo "Welcome!"
echo "This is a script to generate the gnuplot scipt to plot fatbands. You need to provide:"
echo "==========================================="
echo "1. TDOS.dat"
echo "This file is generated from VASPKIT."
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
  else
    file_name="TDOS.dat"
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
# 读取态密度图能量范围
read -p "E_min: " x_min
read -p "E_max: " x_max
# 读取标题
read -p "Title: " temp
if [ "${temp}" != "" ];then
  title=${temp}
else
  title="Total DOS"
fi
##################################################################
##########################   Output part   #######################
##################################################################
echo
echo "Generating *.gnu files ..."
  output="${title}.gnu"
cat >${output}<<- EOF
set encoding iso_8859_1
set terminal  pngcairo truecolor enhanced lw 5.0 font "TimesNewRoman, 44" size 1920, 1680
set output '${title}.png'
set border
set title "${title}" offset 0, -0.8 font "TimesNewRoman, 54"
set style data linespoints
unset ztics
unset key
set view 0,0
set xtics font "TimesNewRoman, 44"
set xtics offset 0, 0.3
set ytics font "TimesNewRoman, 40"
set ylabel font "TimesNewRoman, 48"
set ylabel offset 1.0, 0
set xrange [${x_min}:${x_max}]
set ylabel "DOS"
plot '${file_name}' u (\$1):(\$2+(${de})) w lines lw 1.5 lc "black"
EOF

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
