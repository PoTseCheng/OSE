## OSE Project by Po-Tse Cheng

This project implemented the Generalised Stochastic Simulation Algorithm (GSSA) by [Kenneth L. Judd](https://kenjudd.org/), [Lilia Maliar](https://lmaliar.ws.gc.cuny.edu/), and [Serguei Maliar](https://web.stanford.edu/~maliars/). The original paper: Numerically stable and accurate stochastic simulation approaches for solving dynamic economic models can be found [here](https://onlinelibrary.wiley.com/doi/pdf/10.3982/QE14), while the code of the original GSSA can be found [here](https://web.stanford.edu/~maliars/Files/LowerErrorBounds_ECMA_JMM_2016.zip).

This project translate the entire Matlab code from the original authors to Python code. More importantly, I highlight how GSSA overcomes the two major issues that plague normal stochastic simulation methods: (1) Numerical instability in higher polynomials, and (2) No accuracy improvement in higher polynomials. GSSA conquers the former by using regression methods that can deal with ill-conditioned data, while increase the accuracy by incorporating other integration method. Each method is thoroughly examined and the codes are well-documented. I also include the extension from the author, which solves a multi-country model using GSSA. The results are amazingly accurate. Overall, GSSA is a good choice for solving stochastic economic models.

Please use the following badges to view my notebook:

<a href="https://nbviewer.jupyter.org/github/PoTseCheng/Microeconometrics/blob/master/Final_project.ipynb"
   target="_parent">
   <img align="center"
  src="https://raw.githubusercontent.com/jupyter/design/master/logos/Badges/nbviewer_badge.png"
      width="109" height="20">
</a>
<a href="https://mybinder.org/v2/gh/HumanCapitalAnalysis/microeconometrics-course-project-PoTseCheng.git/master?filepath=Final_project.ipynb"
    target="_parent">
    <img align="center"
       src="https://mybinder.org/badge_logo.svg"
       width="109" height="20">
</a>

Side note: The environment I am using can be found in the environment.yml in this repository. 


---
This Notebook is maintained by Po-Tse Cheng (Viktor Cheng) 
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://github.com/HumanCapitalAnalysis/template-course-project/blob/master/LICENSE)
