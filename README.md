# Dissecting the inactivation mechanism of KcsA with free energy molecular dynamics simulations

**Authors: Sergio Pérez-Conesa and Lucie Delemotte**

------------

<div align="center"><p>
<a href="">
  <img src="https://img.shields.io/badge/python-%2314354C.svg?style=for-the-badge&logo=python&logoColor=white" alt="Python">
</a>
<a href="">
  <img src="https://img.shields.io/badge/Made%20with-Jupyter-orange?style=for-the-badge&logo=Jupyter" alt="Jupyter">
</a>
<a href="">
  <img src="https://img.shields.io/badge/VIM-%2311AB00.svg?style=for-the-badge&logo=vim&logoColor=white" alt="VIM">
</a>
<a href="https://www.linkedin.com/in/sperezconesa/">
  <img src="https://img.shields.io/badge/linkedin-%230077B5.svg?style=for-the-badge&logo=linkedin&logoColor=white" alt="LinkedIn">
</a>
</p>
</div>

<div align="center"><p>
<a href="https://github.com/psf/black">
  <img src="https://img.shields.io/badge/code%20style-black-000000.svg" alt="Code style: black">
</a>
<a href="https://lbesson.mit-license.org/">
  <img src="https://img.shields.io/badge/License-MIT-blue.svg" alt="MIT license">
</a>
<a href="">
  <img src="http://img.shields.io/badge/DOI-XXXXX-B31B1B.svg" alt="DOI:XXXXX">
</a>
<a href="https://twitter.com/intent/follow?screen_name=sperezconesa">
  <img src="https://img.shields.io/twitter/follow/sperezconesa?style=social&logo=twitter" alt="follow on Twitter">
</a>
<a href="https://twitter.com/intent/follow?screen_name=delemottelab">
  <img src="https://img.shields.io/twitter/follow/delemottelab?style=social&logo=twitter" alt="follow on Twitter">
</a>
<a href="https://github.com/sperezconesa/KcsA_string_method_FES">
    <img title="Star on GitHub" src="https://img.shields.io/github/stars/sperezconesa/KcsA_string_method_FES.svg?style=social&label=Star">
</a>
</p>
</div>



------------
![](./reports/final_figures/plots/FES_LB-CHARMM.png)

This project has a description given here by [Sergio Pérez-Conesa](https://www.linkedin.com/in/sperezconesa/). We are  members of the [Delemottelab](https://github.com/delemottelab) led by [prof. Lucie Delemotte](https://www.biophysics.se/index.php/members/lucie-delemotte/). All the explanations can be found in the article and the rest of code and data [here](https://osf.io/snwbc/?view_only=1338fd9e92f941deb7452525c1e9fdfa)

I am happy to connect and discuss this and other projects through [github](https://github.com/sperezconesa), [linkedin](https://www.linkedin.com/in/sperezconesa), [twitter](https://twitter.com/sperezconesa), [email](sperezconesa@gmail.com) etc.
Feel free to suggest ways we could have improved this code.

You can find more updates on the Delemottelab on [twitter](https://twitter.com/delemottelab) and the rest of our [publications](https://scholar.google.es/citations?user=OaHNSvEAAAAJ&hl=en&oi=ao).

If you want to cite this code, please use CITE.bib, thank you!

Published Preprint: Coming soon :wink: []()

Published Article: Coming soon :wink: []()

## Running the code


### Recreate conda environment

To recreate the conda environment used:

```bash
conda env create -f environment.yml
conda activate string_sims
ipython kernel install --user --name=string_sims
pip install -e .
```

Use `environment_exact.yml` for the exact environment.

### Getting additional data files

All the data, including the inference models, simulations etc. can be found in [Open Software Foundation](https://osf.io/snwbc/?view_only=1338fd9e92f941deb7452525c1e9fdfa).

## Project Organization

```text
├── LICENSE
│
├── Makefile           <- Makefile with commands like `make update_data` or `make format`
│
├── README.md          <- The top-level README for developers using this project.
├── data
│   ├── external       <- Data from third party sources.
│   ├── raw            <- Raw data generated or colected for this work
│   ├── interim        <- Intermediate data that has been transformed.
│   └── processed      <- The final, canonical data sets for modeling.
│
├── models             <- MD and experimental models and input files
│
├── notebooks          <- Jupyter notebooks.
│
├── reports            <- Generated analysis as HTML, PDF, LaTeX, etc.
│   ├── tex            <- Latex input files
│   ├── figures        <- Figures of the data.
│   └── final_figures  <- Generated graphics and figures to be used in final report.
│
├── environment.yml    <- The necessary packages to install in conda environment.
│
├── environment_exact.yml   <- The exact package versions used.
│
├── setup.py           <- File to install python code
│
├── src                <- Source code for use in this project.
│   ├── analysis       <- Python code for analysis of data
│   ├── data           <- Python code for handling the data
│   ├── __init__.py    <- Makes src a Python module
│   └── data           <- Scripts to download or generate data
│
├── visualization                <- File to help visualize data and trajectories
│
└──
```

------------

Project based on the [cookiecutter for Molecular Dynamics](https://github.com/sperezconesa/cookiecutter-md). Which is itself based on the [cookiecutter data science project template](https://drivendata.github.io/cookiecutter-data-science/) \#cookiecutterdatascience

------------

## To Do

- [x] Clean-up the repository of unused or WIP files.
- [x] Make a version specific `environment_exact.yml`.
- [ ] Rewrite `README.md`.
- [ ] Save notebooks with images.
- [x] Add github badges.
- [ ] Update github and go public.
- [ ] Update arxiv link.
- [ ] Update `CITE.bib` and doi badge.
- [x] Make repo for big files and link it in `README.md`.
- [ ] Update article link and doi badge.
