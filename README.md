# PFT

This project is a C++ implementation of [Fast and Accurate Partial Fourier Transform for Time Series Data](https://dl.acm.org/doi/10.1145/3447548.3467293) (KDD 2021).
The codes for [FFTW](http://www.fftw.org/index.html), 
[MKL](https://software.intel.com/mkl), 
[Pruned FFT](http://www.fftw.org/pruned.html), and
[Goertzel algorithm](https://github.com/pramasoul/jrand/blob/master/goertzel.c) are also included in `src/`.

## Prerequisites

The implementation requires the following libraries.

- mkl.h
- mkl_dfti.h
- ipp.h
- ipps.h
- fftw3.h

## Datasets

The four datasets used in our paper are available [here](https://drive.google.com/file/d/1ArejxayJdkCTitxhY42iVCCTNpDIc2yd/view?usp=sharing).
They include synthetic random vectors of length of integer power of two.
[Urban sound](https://urbansounddataset.weebly.com/urbansound8k.html) contains 4347 sound recordings in urban environment,
and [Air condition](https://archive.ics.uci.edu/ml/datasets/Appliances+energy+prediction) is composed of 29 time-series vectors of air condition information (e.g., temperature and humidity).
Stock is a new public data we release; it consists of the daily historical stock prices of FANG, 
the four American technology companies Facebook, Amazon, Netflix, and Google.
We collected closing prices adjusted for stock splits, from 2017-01-03 to 2021-01-08.

## Reference

If you use this code, please cite the following paper.
~~~
@inproceedings{park2021fast,
  title={Fast and accurate partial fourier transform for time series data},
  author={Park, Yong-chan and Jang, Jun-Gi and Kang, U},
  booktitle={Proceedings of the 27th ACM SIGKDD Conference on Knowledge Discovery \& Data Mining},
  pages={1309--1318},
  year={2021}
}
~~~
