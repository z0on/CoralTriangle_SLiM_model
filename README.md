
# Coral Triangle SLiM model

### Scripts and datasets for: ###

Matz, M., Treml, E., and Haller, B. Estimating the potential for coral adaptation to global warming across the Indo-West Pacific. *Global Change Biology* 2020 

Input files for SLiM modeling are:
* migration rates: **tri_migration_680.txt**
* population carrying capacities: **tri_popsize_680.txt**
* environment temperature profiles: files in the **environment_profiles** subdir.

To regenerate these files, use **step0_connMap_chooseData_writeEnvs_v4.R**

The rest of R scripts recreate analyses and figires from the paper using raw results that can be downloaded from [here](https://www.dropbox.com/s/012vm9w327ul8hy/RawResults.zip?dl=0) (unzip them into a subdirectory named **RawResults** within the directory cloned from GitHub).

SLiM model codes have extension **.slim**. There are many of them, for settings that cannot be supplied as external arguments. Here is the legend to their filename parts:

* 85: RCP8.5
* 45: RCP4.5
* sin: sinusoidal temperature variation
* rand: sinusoidal + random temperature variation
* novar: no temperature variation
* ljm: low juvenile mortality
* noEvo: no selection response during warming
* noMig: no migration during warming
* noMut: no new mutations during warming

The "main settings" model is **tri85sin_v3.1.1.slim**

To run some or all models, follow this example: 

```
# running 4 replicates of RCP8.5 ("main setting"):
slim -d seed=1 -d popscale=1 -d mutRate=1e-5 -d envirSD=0.5 -d Plasticity=1 -d MutEffect=0.03 -d numQTLs=100 tri85sin_v3.1.1.slim >32tri8a
slim -d seed=2 -d popscale=1 -d mutRate=1e-5 -d envirSD=0.5 -d Plasticity=1 -d MutEffect=0.03 -d numQTLs=100 tri85sin_v3.1.1.slim >32tri8b
slim -d seed=3 -d popscale=1 -d mutRate=1e-5 -d envirSD=0.5 -d Plasticity=1 -d MutEffect=0.03 -d numQTLs=100 tri85sin_v3.1.1.slim >32tri8c
slim -d seed=4 -d popscale=1 -d mutRate=1e-5 -d envirSD=0.5 -d Plasticity=1 -d MutEffect=0.03 -d numQTLs=100 tri85sin_v3.1.1.slim >32tri8d

# running 4 replicates of RCP4.5 :
slim -d seed=1 -d popscale=1 -d mutRate=1e-5 -d envirSD=0.5 -d Plasticity=1 -d MutEffect=0.03 -d numQTLs=100 tri45sin_v3.1.1.slim >32tri4a
slim -d seed=2 -d popscale=1 -d mutRate=1e-5 -d envirSD=0.5 -d Plasticity=1 -d MutEffect=0.03 -d numQTLs=100 tri45sin_v3.1.1.slim >32tri4b
slim -d seed=3 -d popscale=1 -d mutRate=1e-5 -d envirSD=0.5 -d Plasticity=1 -d MutEffect=0.03 -d numQTLs=100 tri45sin_v3.1.1.slim >32tri4c
slim -d seed=4 -d popscale=1 -d mutRate=1e-5 -d envirSD=0.5 -d Plasticity=1 -d MutEffect=0.03 -d numQTLs=100 tri45sin_v3.1.1.slim >32tri4d

# removing unnecessary lines from output
for F in 32tri*[abcd] ; do 
cat $F | grep -v adults | grep -v empty | grep -v extinct >$F.clean;
done

# place resulting *.clean files into ./RawResults subdirectory
```

