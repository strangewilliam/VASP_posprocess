POTCAR 生成脚本

```shell
#!/bin/bash
# This is potmake script.

hom="/fs12/home/zhj_xujq"
file_name=`find -name  POTCAR*`
if [ "${#file_name[*]}" != "0" ];
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
  read method
  echo
  if [[ ${method} == "LDA" || ${method} == "PBE" ]];
  then
    judge="1"
  else
    echo "Wrong input!"
  
  fi
done
hom="${hom}/POTCAR/${method}.54/"
for ele in ${elements}
do
# judge whether the POTCAR needed exist first
  echo "The types of the POTCAR files you may want:"
  ls -d ${hom}/${ele}*
  echo
  read -p "Please type the file name you want:"  file_name 
  if [ ! -n "${file_name}" ];then
    file_name=${ele}
    echo ${file_name}
  fi
  echo
  while [ ! -e ${hom}/${file_name} ]
  do
    echo "No such file! Please input again!"
    echo
    echo "The types of the POTCAR files you may want:"
    ls -d ${hom}/${ele}*
    read -p "Please type the file name you want:"  file_name 
    if [ ! -n "${file_name}" ];then
      file_name=${ele}
      echo ${file_name}
    fi
    echo
  done
  cat ${hom}/${file_name}/POTCAR >> ./POTCAR_${method}
done

## Find the maximum of the ENMAXs
grep 'ENMAX' POTCAR_${method}>temp
row=`awk 'END{print NR}' temp`
m=0
for r in $( seq ${row} )
do
  t=`awk NR!=${r}'{next}{printf("%.3f",$3)}' temp`
  judge=`echo "${t}>${m}"|bc`
  if [ "${judge}" = "1" ];then
    m=${t}
  fi	
done
rm temp

echo -e "Do you want to rename the POTCAR_${method} file as POTCAR?\nType [y] to rename the POTCAR file or any other key to keep the name or [Ctrl+c] to quit."
read judge
echo
if [ "${judge}" = "y" ];then
  mv POTCAR_${method} POTCAR
  file_name="POTCAR"
else
  file_name="POTCAR_${method}"
fi
num_ele=`grep -i vr ${file_name}|awk 'END{print NR}'`
num=`echo ${elements}|awk 'END{print NF}'`
if [ "${num_ele}" != "${num}" ];then
  echo "The number of elements contained in the ${file_name} file is not correct. May you good luck!"
  rm POTCAR*
else
  echo "ENMAX among elements " ${elements} " is " ${m}
  echo "Successfully finished!"
fi
exit 0
```

