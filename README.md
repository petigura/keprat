# Kepler Radius Ratio

## Cookbook

### Set appropriate conda environment

Sadly, I never got around to coming up with a unique conda environment
for this code.

- Notebooks 1-4 are run in my base python2 environment, which is
  dangerous because I could update it

- Notebooks 5-6 are run in py37, because they need python3 to run.

## Generate plots for paper

run the numberered ipython notebooks in `notebooks/`

## Notes

### How to read in Kepler project chains

```
df = pd.read_hdf('data/kepler_project_chains.hdf','chains-dr25-K00063.01')
```

and column definitions are here

data/kepler-project-mcmc-column-definitions.txt

### Final sample of planets

```
data/cksgaia-planets-tau.csv # 901 planet sample
data/column-definitions.txt # info on the columns
notebooks/1_Fig-Compare-Van-Eylen.ipynb # shows how to plot to make cleaned CKS sample
```
