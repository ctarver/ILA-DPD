# GMP DPD Using Indirect Learning Architecture Matlab Library
[![DOI](https://zenodo.org/badge/142376314.svg)](https://zenodo.org/badge/latestdoi/142376314)

## How to Cite
Please cite the repo if you use it in your projects. An example bibtex entry is below.

```
@misc{TarverILADPD,
  author       = {Tarver, Chance},
  title        = {GMP DPD Using Indirect Learning Architecture Matlab Library},
  month        = apr,
  year         = 2019,
  doi          = {10.5281/zenodo.2648687},
}
```


## Introduction
Power amplifiers (PAs) are nonlinear devices. These nonlinearities contribute to distortions such as spectral regrowth around a carrier. There are limits placed by the 3GPP and FCC on this adjacent carrier leakage and other nonlinear artifacts.

To counter this, we can predistort with an inverse of the PA's nonlinearities. A standard way of doing that is via an indirect learning architecture. A challenge with designing a predistorter is that we don't know what the predistorter output should be, so we can directly use least squares to solve for the predisorter design. An indirect learning architecture allows us to circumvent this.

![ILA](http://zone.ni.com/images/reference/en-XX/help/374264E-01/dpd2.png "Indirect Learning Architecture")

## How to install this: 
#### Option 1: Add ILA_DPD.m to the path:
Download this repo. Put the ILA_DPD.m in the folder of your project that needs DPD or in any folder that is in the Matlab path. Then proceed to use the class to train and use a predistorter.

#### Option 2: Add as a submodule in your git project:
If you already are using git with your project, you can use this as a submodule. In the main directory of your project, run
```
git submodule add https://github.com/ctarver/ILA-DPD.git
git submodule update --init --recursive
```
This repo should show up in your project. Then just commit the new submodule to your project like you would commit anything. 
To use the class, you still need to add the ILA_DPD.m to your path. This can be done via `addpath(genpath('ILA-DPD'))` if your main script is in a directory above the submodule.

## How to use the DPD class:
Initialize with:
`dpd = ILA_DPD(params);`

The 'params' stuct needs to contain the following:
 - `params.order`: Nonlinear order of the DPD model. Must be an odd integer.
 - `params.memory_depth`: Memory depth used on each branch of the parllel hammerstein DPD model. Can be any positive integer.
 - `params.lag_depth`: Memory lag/lead depth used on each cross term branch of the GMP DPD model. Can be any positive integer. If 0, then this becomes a parallel Hammerstein memory polynomial.
 - `params.nIterations`: Number of times to go through the indirect learning iterations. Can be any positive integer. Usually doesn't need to be higher than 4.
  
The DPD model can be trained by `dpd.perform_learning(tx_data, board);`. The `tx_data` is a column vector of IQ data to send through the PA for generating a model. The `board` is a PA class that needs to have a `transmit` method. This works with my [PA model](https://github.com/ctarver/Power-Amplifier-Model) or [WARPLab Wrapper](https://github.com/ctarver/WARPLab-Matlab-Wrapper).

Once the DPD is trained, you can predistort any samples with the model using the `predistort` method. Just do `dpd.predistort(tx_data);`

Check out the example.m to see more. 

## Result: 
Running this code on RFWebLab, I obtained the following power spectral density plot. Here, I am using a 20 MHz LTE like signal. With DPD, we can see less spectral regrowth around the carrier. 
Here, ***P*** represents the nonlinearity order, ***M*** is the memory depth, and ***L*** is the lag/lead amount. 
![PSD](https://raw.githubusercontent.com/ctarver/ILA-DPD/master/psd_example.png?token=ACLnMTVWU6jnNqXKfcndnWRs5eeq5Ph8ks5bZG90wA%3D%3D "RFWebLab PSD")

## Theory of Operation:
The DPD is a nonlinear function that approximates an inverse of the PA's nonlinearities so that the sum output is linear. The DPD in principle can be modeled in any fashion. I have chosen to use a Parallel Hammerstein model. This consists of multiple nonlinear branches where each is followed by an FIR filter to model memory effects. The parallel Hammerstein model is widely chosen for its balance between modeling performance and complexity. 

One strong advantage of this model is that the FIR filter and hence any coefficients of the model are after their corresponding nonlinearities. This means the output of the model can be expressed as a linear system, *y = X b*. Here, X is a matrix where each column is a different nonlinear branch of the form *x|x|^{i-1}* where *i* is the order which can only be odd, and *x* is the input signal. If we have input/output samples, this can be a linear least-squares regression problem since the model coefficients are linear with regards to the input. Here we want to *minimize_b || y - X b ||* for some experimental *x* and *y.* This minimizes the sum of squared residuals and gives us the best fit. 

The DPD can't easily be modeled in this method directly because we don't know the desired predistorter output. We just know the desired PA output, actual PA output, and the original input signal. The indirect learning architecture allows us to circumvent this. When fully converged, the PA output should be linear, and so the input to the pre and post distorters would be equivalent, their outputs would be the same, and the error signal would be zero. When training, we use the postdistorter.  We want the output of the postdistorter to be equal to the output of the predistorter so that there is no error. We can run this to train and get some PA output signal. Then for the postdistorter, we have input sample (the PA ouput) and a desired postdistorter output (the predistorter output). We can start with some DPD coefficients (such as a pure linear DPD) then perform a LS fit to find the best coefficients to fit the postdistorter.  We copy this to the predistorter and repeat for a few iterations.

## References:
For Generalized memory polynomials:
```
D. R. Morgan, Z. Ma, J. Kim, M. G. Zierdt and J. Pastalan, "A Generalized Memory Polynomial Model for Digital Predistortion of RF Power Amplifiers," in IEEE Transactions on Signal Processing, vol. 54, no. 10, pp. 3852-3860, Oct. 2006.
doi: 10.1109/TSP.2006.879264
```

For the conjugate branch:
```
L. Anttila, P. Handel and M. Valkama, "Joint Mitigation of Power Amplifier and I/Q Modulator Impairments in Broadband Direct-Conversion Transmitters," in IEEE Transactions on Microwave Theory and Techniques, vol. 58, no. 4, pp. 730-739, April 2010.
doi: 10.1109/TMTT.2010.2041579
keywords: {modulators;power amplifiers;radio transmitters;telecommunication channels;broadband direct-conversion transmitters;frequency-dependent power amplifier;I/Q modulator impairments;direct-conversion radio transmitters;extended parallel Hammerstein structure;parameter estimation stage;indirect learning architecture;adjacent channel power ratio;Broadband amplifiers;Power amplifiers;Radio transmitters;Digital modulation;Predistortion;Local oscillators;Wideband;Nonlinear distortion;Radiofrequency amplifiers;Frequency estimation;Digital predistortion (PD);direct-conversion radio;in-phase and quadrature (I/Q) imbalance;I/Q modulator;local oscillator (LO) leakage;mirror-frequency interference (MFI);power amplifier (PA);spectral regrowth},
URL: http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=5431085&isnumber=5446455
```

## History

### v1.1
Sept 27th, 2019

* Now works for GMP and includes a lag/lead crossterm.

### v1.0
April 2019

* First official tagged release with a DOI. 
* Works for standard memory polynomials.
* Supports conjugate branch processing for IQ imbalance
