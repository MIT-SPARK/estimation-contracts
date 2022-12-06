# README

# Overview

This is code to reproduce the experiments in the paper: L. Carlone, [**Estimation Contracts for Outlier-Robust 
Geometric Perception**], [arXiv: 2208.10521](https://arxiv.org/pdf/2208.10521.pdf), 2022].

# Getting started

The code has been implemented and tested in Matlab R2020a, using a Macbook Pro with maxOS Monterey. We expect to run on any other operating system supporting the dependencies below.

To run the code, please install the following dependencies. Note: it is important to install all the dependencies in the same folder (e.g., all subfolders of a ``code'' folder) 
- [CVX with MOSEK](http://cvxr.com/cvx/doc/mosek.html) (tested with CVX Version 2.2)
- [STRIDE and the Certifiably Robust Perception repository](https://github.com/MIT-SPARK/CertifiablyRobustPerception/) (tested with latest commit on December 10, 2022)
- [SOSTOOLS](https://github.com/oxfordcontrol/SOSTOOLS) (tested with Version 3.01)

Clone this repo in the same folder as the other dependencies. 
After installing the dependencies and cloning this repo, the folder where you installed the dependencies should look like this:
-

## Reference

If you found the paper or code useful, please cite:

```bibtex
@article{Carlone22arxiv-estimationContracts,
  author = {L. Carlone},
  title = {Estimation Contracts for Outlier-Robust Geometric Perception},
  journal = {arXiv preprint arXiv: 2208.10521},
  note = {\linkToPdf{https://arxiv.org/pdf/2208.10521.pdf}},
  pdf = {https://arxiv.org/pdf/2208.10521.pdf},
  Year = {2022}
}
```

## Acknowledgments

This work was partially funded by the NSF CAREER award [Certifiable Perception for Autonomous Cyber-Physical
Systems](https://nsf.gov/awardsearch/showAward?AWD_ID=2044973) and by ARL [DCIST](https://www.dcist.org/) CRA W911NF-17- 2-0181.

## License

[BSD License](LICENSE.BSD)



