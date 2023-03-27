# Toolbox for fast, variational 1D ERP alignment (VERPA)

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

If you use this code for other work, please additionally cite
  
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

## Toolbox for nonlinear, anisotropic diffusion on the GPU
Fast and fully vectorized MATLAB implementation of non-linear, anisotropic diffusion with gpuarray and RGB image support. 
Implemented diffusion types are edge enhancing and coherence enhancing diffusion (Weickert, 1999) as well as Perona-Malik non-linear isotropic diffusion with different diffusivities. 

Download
```
$ git clone https://github.com/phflot/diffusion_toolbox.git
```

Documentation and Usage
The file ```demo.m``` demonstrates the usage of the functions. 

## ASPR ERP Image Denoising
if you use ASPR image denoising for your work, please cite 

Reference
```
@article{kohl19,
  author  = {Kohl, M. and Schebsdat, E. and Schneider, E. and Strauss, D. J.},
  journal = {Conference proceedings: 2019 9th International IEEE/EMBS Conference on Neural Engineering (NER)},
  title   = {Denoising of Single-Trial Event-Related Potentials by Shrinkage and Phase Regularization of Analytic Wavelet Filterbank Coefficients},
  year    = {2019},
  pages   = {251-54},
  doi     = {10.1109/NER.2019.8717148},
}
```

This repository provides a self-contained MATLAB implementation of an algorithm for the fast denoising of single-trial event-related 

Download
```
$ git clone https://gitlab.com/manuelchristophkohl/aspr-erp-image-denoising.git
```

If you use the BM3D functions in this repository for your work, please cite

## BM3D demo software for image/video restoration and enhancement  
Public release v1.9 (26 August 2011) 

Download
```
https://webpages.tuni.fi/foi/GCF-BM3D
```

Reference
```
@article{Dabov2007,
  author  = {Dabov, K. and Foi, A. and Katkovnik, V. and Egiazarian, K.},
  journal = {IEEE Trans Image Process},
  title   = {Image Denoising by Sparse 3-D Transform-Domain Collaborative Filtering},
  year    = {2007},
  pages   = {2080-95},
  volume  = {16},
  doi     = {10.1109/TIP.2007.901238},
}
```

## Anisotropic phase-map denoising
If you use the Anisotropic phase-map denoising function in this repository for your work, please cite

Reference
```
@article{Villa2010,
  author   = {Villa, J. and Rodriguez-Vera,R. and Quiroga,J. A. and de la Rosa, I. and González,E.},
  journal  = {Opt Lasers Eng},
  title    = {Anisotropic phase-map denoising using a regularized cost-function with complex-valued Markov-random-fields},
  year     = {2010},
  pages    = {650-56},
  volume   = {48},
  doi      = {https://doi.org/10.1016/j.optlaseng.2010.02.002},
  keywords = {Regularization, Phase-map filtering},
  url      = {http://www.sciencedirect.com/science/article/pii/S0143816610000242},
}
```

## Woody and thornton average by
Ikaro Silva (2023). Woody Average (Average with alignment)
(https://www.mathworks.com/matlabcentral/fileexchange/12459-woody-average-average-with-alignment), MATLAB Central File Exchange
