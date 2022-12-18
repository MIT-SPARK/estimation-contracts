# README

# Overview

This is the code to reproduce the experiments in the paper: 

**L. Carlone, Estimation Contracts for Outlier-Robust Geometric Perception, [arXiv: 2208.10521](https://arxiv.org/pdf/2208.10521.pdf), 2022.**

# Getting started

The code has been implemented and tested in Matlab R2020a, using a Macbook Pro with maxOS Monterey. We expect it to run on any other operating system supporting the dependencies below.

To run the code, please install the following dependencies. Note: it is important to install all the dependencies in the same folder (e.g., all subfolders of a ``code`` folder):
- [CVX with MOSEK](http://cvxr.com/cvx/doc/mosek.html) (tested with CVX Version 2.2)
- [STRIDE and the Certifiably Robust Perception repository](https://github.com/MIT-SPARK/CertifiablyRobustPerception/) (tested with latest commit, December 10, 2022)
- [SOSTOOLS](https://github.com/oxfordcontrol/SOSTOOLS) (tested with Version 3.01/master branch, but also runs with SOSTOOLS400 branch, which is typically faster)

Clone our repo in the same folder as the other dependencies. 
After installing the dependencies and cloning this repo, the folder where you installed the dependencies in should at least have the following sub-folders:
- CertifiablyRobustPerception
- cvx
- estimation-contracts
- SOSTOOLS

# Running the code

You can replicate the experiments in the paper by running each of the following experiments:
* **experiment1_aposteriori_bounds_synthetic**: solves synthetic rotation search problems using [QUASAR](https://arxiv.org/pdf/1905.12536.pdf) and compares the estimation errors with the a posteriori bound developed in our paper (Fig. 4 in the paper) (expected runtime: ~15 hours)
* **experiment2_check_hypercontractivity_synthetic**: check certifiable hypercontractivity in synthetic rotation search problems (Fig. 5 in the paper) (expected runtime: ~1.5 hours)
* **experiment3_check_hypercontractivity_real**: check certifiable hypercontractivity in real rotation search problems arising in panorama stitching (Table 1 in the paper). To run the real tests, please download the data from and put them in the **data** subfolder: [link to panorama stitching data](https://drive.google.com/drive/folders/1CppsDdU98PgG939aV0ZaaBcVYRLrgI9O?usp=sharing) 
(expected runtime: ~0.5 hours)
* **experiment4_hypercontractivity_bounds**: visualizes the bounds in Theorem 11 (Fig. 6 in the paper)
* **experiment5_visualize_anticoncentration**: visualizes the functions involved in the definition of certifiable anti-concentration (Fig. 7 in the paper)
* **experiment6_check_anticoncentration_synthetic_vs_n**: check certifiable anti-concentration in synthetic rotation search problems for different number of measurements (Fig. 8(b) in the paper) (expected runtime: ~55 hours)
* **experiment7_check_anticoncentration_synthetic_vs_eta**: check certifiable anti-concentration in synthetic rotation search problems for different values of the parameter $\eta$, which obtains:
  * Fig. 8(a) when the parameter ``isModified = 0`` (expected runtime: ~40 hours)
  * Fig. 8(c) when the parameter ``isModified = 0`` (expected runtime: ~40 hours)
* **experiment8_slide_synthetic**: evaluates SLIDE for list decodable regression in synthetic rotation search problems, which obtains:
  * Fig. 9 when the parameter ``isAdversarial = 0`` and ``recoverAllHypotheses = 0`` (expected runtime: ~6 hours)
  * Fig. 10 when the parameter ``isAdversarial = 1`` and ``recoverAllHypotheses = 0`` (expected runtime: ~6 hours)
  * Fig. 11 when the parameter ``isAdversarial = 1`` and ``recoverAllHypotheses = 1`` (expected runtime: ~6 hours)
  * Fig. 12 when the parameter ``isAdversarial = 1`` and ``recoverAllHypotheses = 1`` (for some random instances of the problem, also stored in SLIDE_example1.mat and SLIDE_example2.mat)
* **experiment9_check_anticoncentration_real**: check certifiable anti-concentration in real rotation search problems arising in panorama stitching (not reported in the paper). To run the real tests, please download the data from and put them in the **data** subfolder: [link to panorama stitching data](https://drive.google.com/drive/folders/1CppsDdU98PgG939aV0ZaaBcVYRLrgI9O?usp=sharing)

The results we obtained by running the examples above are also stored in the **results** folder within this repo.
Most results are obtained using SOSTOOLS v. 3.01, but we re-ran the experiments where we measured the runtime using the SOSTOOLS400 branch in SOSTOOLS, since the latter is typically faster. We store both sets of results and add a label *SOSTOOLS400* to the latter set of results.

# Reference

If you found the paper or code useful, please cite:

```bibtex
@article{Carlone22arxiv-estimationContracts,
  author = {L. Carlone},
  title = {Estimation Contracts for Outlier-Robust Geometric Perception},
  journal = {arXiv preprint arXiv: 2208.10521},
  pdf = {https://arxiv.org/pdf/2208.10521.pdf},
  Year = {2022}
}
```

# Acknowledgments

This work was partially funded by the NSF CAREER award [Certifiable Perception for Autonomous Cyber-Physical
Systems](https://nsf.gov/awardsearch/showAward?AWD_ID=2044973) and by ARL [DCIST](https://www.dcist.org/) CRA W911NF-17- 2-0181.

# License

[BSD License](LICENSE.BSD)



