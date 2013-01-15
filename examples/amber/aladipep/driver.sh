for i in `seq 1 21`
do 
echo $i 

less pathdata.temp | sed s/TOSET/$i/g > pathdata 

sh run.sh 

done 

