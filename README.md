# KEMAR HRTFs with low frequency correction
This repository holds HRTFs of a KEMAR 45BA dummy head with large ears. The HRTFs are a modified version from [previous measurements](https://github.com/sfstoolbox/data/blob/master/HRTFs/QU_KEMAR_anechoic_3m.sofa) and have been created by the Matlab script [Correct_low_frequencies_of_HRTFs.m](http://github.com/spatialaudio/hptf-compensation-filters/blob/master/Correct_low_frequencies_of_HRTFs.m).

## Method for low frequency correction of HRTFs
Based on
* Xie, B. (2009): On the low frequency characteristics of head-related transfer functions. Chinese J. Acoust. 28(2), pp. 1-13

the magnitude response of the HRTFs has been set to a constant value and the phase has been extrapolated linearly for low frequencies. As the group delay at low frequencies is now low enough, the HRTFs can be truncated to 512 samples.
Alternative methods for low frequency correction of HTRFs can be found in
* Algazi, V. R., Duda, R. O., Duraiswami, R., Gumerov, N. A. and Tang, Z. (2002): Approximating the head-related transfer function using simple geometric models of the head and torso. J. Acoust. Soc. Am. 112(5), pp. 2053-2964
* Gumerov, N. A., O'Donovan, A. E., Duraiswami, R. and Zotkin, D. N. (2010): Computation of the head-related transfer function via the fast multipole accelerated boundary element method and its spherical harmonic representation. J. Acoust. Soc. Am. 127(1), pp. 370-386
* Bernschütz, B. (2013): A Spherical Far Field HRIR/HRTF Compilation of the Neumann KU 100. Proc. of the 39th German Annual Conf. on Acoust. (DAGA)

## Licenses
The Matlab code by Vera Erbes is licensed under the MIT license. Confer header of [Correct_low_frequencies_of_HRTFs.m](http://github.com/spatialaudio/hptf-compensation-filters/blob/master/Correct_low_frequencies_of_HRTFs.m) for more information on this license.
