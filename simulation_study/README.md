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
 └── binary_estimation: Folder containing simulation code and results for assessing estimation of the overall treatment effect in simulated datasets with a binary outcome
    ├── README.md
    ├── bash.sh
    ├── main.R
    └── results 
└── cont_estimation: Folder containing simulation code and results for assessing estimation of the overall treatment effect in simulated datasets with a continuous outcome
    ├── README.md
    ├── bash.sh
    ├── main.R
    └── results
├── helpers: helper functions for binary_estimation and cont_estimation
│   ├── helper.R
    ├── README.md
├── parse_simulation_results.R: Code to parse simulation results and make tables shown in paper
└── tests: Folder containing simulation code and results for assessing test Type I error and power (permutation, LR tests)
    ├── README.md
    ├── bash.sh
    ├── main.R
    └── results


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
