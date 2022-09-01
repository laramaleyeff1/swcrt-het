<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#about-the-project">About The Project</a>
    </li>
    <li><a href="#contact">Contact</a></li>
    <li><a href="#acknowledgments">Acknowledgments</a></li>
  </ol>
</details>



<!-- ABOUT  -->
## About The Project

This folder contains the simulation code and results for “Assessing exposure-time treatment effect heterogeneity in stepped wedge cluster randomized trials.” We also provide the code to create the simulation tables and figures shown in the paper.

    .
    ├── README.md
    ├── binary_expt_estimation: Folder containing simulation code and results for assessing estimation of the exposure-time specific 
    effects in simulated datasets with a binary outcome
        ├── results
        ├── results_no_boots
        ├── README.md
        ├── bash.sh 
        ├── bin_boots_10.Rda
        ├── main.R                  
        └── param_no_boots.Rda
    ├── binary_overall_estimation: Folder containing simulation code and results for assessing estimation of the average treatment 
    effect in simulated datasets with a binary outcome
        ├── results
        ├── README.md         
        ├── bash.sh                     
        ├── bin_all_less.Rda
        └── main.R  
    ├── cont_overall_estimation: Folder containing simulation code and results for assessing estimation of the average treatment 
    effect in simulated datasets with a continuous outcome
        ├── results
        ├── README.md 
        ├── bash.sh 
        ├── cont_all.Rda
        └── main.R
    ├── helpers: helper functions for binary_expt_estimation, binary_overall_estimation, and cont_estimation
        ├── README.md
        └── helper.R
    ├── tests: Folder containing simulation code and results for assessing test Type I error and power (permutation, LR tests)
        ├── results
        ├── README.md 
        ├── bash.sh 
        ├── main.R
        └── power.Rda
    └── parse_simulation_results.R: Code to parse simulation results and make tables shown in paper



<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- CONTACT -->
## Contact

Lara Maleyeff - lmaleyeff@g.harvard.edu

Project Link: [https://github.com/laramaleyeff1/swcrt-het](https://github.com/laramaleyeff1/swcrt-het)

<p align="right">(<a href="#readme-top">back to top</a>)</p>


<!-- ACKNOWLEDGMENTS -->
## Acknowledgments

I would like to thank my co-authors Rui Wang, Fan Li and Sebastien Haneuse.

<p align="right">(<a href="#readme-top">back to top</a>)</p>
