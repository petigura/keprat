## Kepler Radius Ratio

### Conda environment

Work in petigura's base conda enviroment (dangerous b/c should update it)


### Reading in Kepler project chains

```
df = pd.read_hdf('data/kepler_project_chains.hdf','chains-dr25-K00063.01')
```

Column defintions are in 

```
data/kepler-project-mcmc-column-definitions.txt
```

###

Create `data/cksgaia-planets.csv` for Vincent to read in with 


```
bin/run_keprat.py cksgaia-planets
```

## To generate plots for papers run the numberered ipython notebooks



