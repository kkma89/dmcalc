DMCalc
======

`DMCalc` estimates the Dispersion Measure (DM) of wide-band pulsar data in [psrfits][psrfits] format. It uses the tools available with [PSRCHIVE][psrchive] python interface to get ToAs and [TEMPO2][tempo2] for DM fitting. A median absolute deviation (MAD) based ToA rejection algorithm following Tiburzi et al. (2019) is also implemented in the code to remove large outlier ToAs. 

## Dependencies
DMCalc is a python

## Using DMCalc


[psrfits]: https://www.atnf.csiro.au/research/pulsar/psrfits_definition/Psrfits.html
[psrchive]: http://psrchive.sourceforge.net/
[tempo2]: https://bitbucket.org/psrsoft/tempo2/src/master/
