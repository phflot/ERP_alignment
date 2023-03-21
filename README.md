# Toolbox for fast, variational 1D ERP alignment

The code is written in Matlab and C++

## Download

Download the repository via
```
$ git clone https://github.com/phflot/ERP_alignment.git
```


## Documentation and Usage

The Code can be compiled and installed with the ```set_path.m``` file. 

## Citation

If you use this work in the context of EP and ERP analysis, please cite

> D. Thinnes, P. Flotho, F. I. Corona-Strauss, D. J. Strauss and J. F. Vibell, “Compensation of P300 Latency Jitter using fast variational 1D Displacement Estimation” (in preparation), 2023. 

BibTeX entry
```
@article{thinn23,
    author = {Thinnes, D. and Flotho, P. and Corona-Strauss, F. I. and Strauss, D. J. and Vibell, J. F.},
    title = {Compensation of P300 Latency Jitter using fast variational 1D Displacement Estimation},
    journal = {(in preparation)},
    year = {2023}
}
```

If you use this code for your other work, please additionally cite
  
> Flotho, P., Thinnes, D., Kuhn, B., Roome, C. J., Vibell, J. F., & Strauss, D. J. (2021). Fast variational alignment of non-flat 1D displacements for applications in neuroimaging. Journal of Neuroscience Methods, 353, 109076.

BibTeX entry
```
@article{floth21,
  title={Fast variational alignment of non-flat 1D displacements for applications in neuroimaging},
  author={Flotho, Philipp and Thinnes, David and Kuhn, Bernd and Roome, Christopher J and Vibell, Jonas F and Strauss, Daniel J},
  journal = {Journal of Neuroscience Methods},
  volume = {353},
  pages = {109076},
  year = {2021},
  issn = {0165-0270},
  doi = {https://doi.org/10.1016/j.jneumeth.2021.109076},
  url = {https://www.sciencedirect.com/science/article/pii/S016502702100011X},
}
```



## Denoising algorithms
If you use nonlinear, anisotropic diffusion for your work, please cite

# Toolbox for nonlinear, anisotropic diffusion on the GPU
Fast and fully vectorized MATLAB implementation of non-linear, anisotropic diffusion with gpuarray and RGB image support. 
Implemented diffusion types are edge enhancing and coherence enhancing diffusion (Weickert, 1999) as well as Perona-Malik non-linear isotropic diffusion with different diffusivities. 

Download
Download the repository via
```
$ git clone https://github.com/phflot/diffusion_toolbox.git

Documentation and Usage
The file ```demo.m``` demonstrates the usage of the functions. 

```
## ASPR ERP Image Denoising
if you use ASPR image denoising for your work, please cite 

This repository provides a self-contained MATLAB implementation of an algorithm for the fast denoising of single-trial event-related 

Download
$ git clone https://gitlab.com/manuelchristophkohl/aspr-erp-image-denoising.git

```
If you use the BM3D functions in this repository for your work, please additionally cite
```
# BM3D demo software for image/video restoration and enhancement  
Public release v1.9 (26 August 2011) 

Download
https://webpages.tuni.fi/foi/GCF-BM3D

Reference
K. Dabov, A. Foi, V. Katkovnik, and K. Egiazarian, "Image 
denoising by sparse 3D transform-domain collaborative filtering," 
IEEE Trans. Image Process., vol. 16, no. 8, August 2007.
