### 0: Initial Set-Up 

Cloning the repository:
```
git clone https://gitlab.cern.ch/shbansal/hplusWh_qqbb.git
```
You would have to change the names of the paths in the RunJob.sh, submitHplusjobs.sh, changing the paths where you store the ntuples etc in diffferent files.

### 1: Running the Code
I usually first make this code in an interactive job on BAF, i.e do
```
condor_submit -interactive CentOS7_interactive.jdl
```
and then switching to the directory where you have stored this code
```
setupATLAS
lsetup "root 6.18.04-x86_64-centos7-gcc8-opt"
make
```
For submitting the jobs, switch to a normal working area within your code directory (i.e in a different terminal/shell --not interactive--):
```
source submitHplusjobs.sh
```