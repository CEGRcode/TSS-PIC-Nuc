module load anaconda3
source activate bioinfo
# Script for F2a and E5

### CHANGE ME
WRK=/Path/to/Title/

###SCRIPT
SCRIPTMANAGER=$WRK/bin/ScriptManager-v0.15.jar
COMPOSITEFILTER=$WRK/bin/sum_Col_CDT_filter.pl
## determin output
[ -d logs ] || mkdir logs
[ -d $WRK/Library/F3f ] || mkdir -p $WRK/Library/F3f
[ -d $WRK/Library/E10 ] || mkdir -p $WRK/Library/E10

cd $WRK/Library/F3f

cd $WRK/Library/E10

