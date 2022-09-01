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

This folder contains the simulation code and results for all continuous outcome simulations in “Assessing exposure-time treatment effect heterogeneity in stepped wedge cluster randomized trials.” We also provide the code to create the simulation tables and figures shown in the paper.

        .
        ├── README.md
        ├── bash.sh                     # Bash script to run from computing cluster, to run: sbatch bash.sh B n_per t_max each sd_expt         
        ├── main.R                      # The main code of the simulation study. Generates a stepped wedge CRT with a variety of calendar and exposure time
                                        # patterns and fits Models 1-5. Outputs results to "results" folder. Assumes that first group of clusters crossover
                                        # in the second time period, with an equal number crossing over in each time period t, for t>=2.

                                        # The scenarios used in the code are as follows:
                                          # Scenario 1 = no exposure time treatment effect heterogeneity
                                          # Scenario 10 = normally distributed exposure time treatment effect heterogeneity
                                          # Scenario 2 = linearly increasing exposure time treatment effect heterogeneity
                                          # Scenario 4 = delayed exposure time treatment effect heterogeneity
        ├── bin_boots_10.Rda            # Results for Scenario 10 only, with standard errors computed using within-cluster bootstrap
        ├── param_no_boots.Rda          # Results for all scenarios - please ignore standard error estimates (depreciated)
        ├── results                     # Raw csv results for Scenario 10 only, with standard errors computed using within-cluster bootstrap
        └── results_no_boots            # Raw csv results for all scenarios - please ignore standard error estimates (depreciated)


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
