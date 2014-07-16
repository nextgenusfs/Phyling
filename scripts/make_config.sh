module load EMBOSS
module load hmmer/3.0
module load cdbfasta
module load fastx_toolkit
module load wise
module load FastTree
module load muscle
module load trimal
module load phrap

TEMPLATE=lib/apps.conf.template
CONF=lib/apps.conf
rm -f $CONF
for r in `sort $TEMPLATE`
do
 VAR=`echo $r | awk -F= '{print $1}'`
 VAL=`echo $r | awk -F= '{print $2}'`
 NEWVAL=`which $VAL`
 echo "$VAR=$NEWVAL" >> $CONF
done
