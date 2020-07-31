# Kepler Radius Ratio

## Cookbook

### Set appropriate conda environment

Sadly, I never got around to coming up with a unique conda environment
for this code.

- Notebooks 1-4 are run in my base python2 environment, which is
  dangerous because I could update it

- Notebooks 5-6 are run in py37, because they need python3 to run.

### Generate plots for paper

run the following notebooks in `in notebooks`

```bash
1_Fig-Compare-Van-Eylen.ipynb
2_Tau-Cut-Simulations.ipynb
3_number-of-short-cadence.ipynb
4_Violin-Plots.ipynb
5_fits-to-simulated-data.pdf.ipynb
6_fits-to-simulated-data_per=115.ipynb
```

## Notes

### How to read in Kepler project chains

```python
df = pd.read_hdf('data/kepler_project_chains.hdf','chains-dr25-K00063.01')
```

and column definitions are here

```bash
data/kepler-project-mcmc-column-definitions.txt
```

### Final sample of planets

```bash
data/cksgaia-planets-tau.csv # 901 planet sample
data/column-definitions.txt # info on the columns
notebooks/1_Fig-Compare-Van-Eylen.ipynb # shows how to plot to make cleaned CKS sample
```
