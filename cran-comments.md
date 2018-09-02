## Test environments
* local ubuntu 18 install, R 3.4.4
* local MS Windows, R 3.3.0 with Rtools
* ubuntu 14.04 (on travis-ci), R 3.4.4, 3.5.0, 3.6.0
* OSX 10.11 (on travis-ci), R 3.4.4, 3.5.0, 3.6.0
* MS Windows x86_64-w64-mingw32/x64 (on appveyor), R 3.5.1

## R CMD check results
There were no ERRORs or WARNINGs. 

There were 2 NOTEs:
* checking CRAN incoming feasibility ... NOTE
Maintainer: ‘Vojtech Spiwok <spiwokv@vscht.cz>’

New submission

* checking examples ... NOTE
Examples with CPU or elapsed time > 5s
              user system elapsed
lines.nebpath 5.04  0.004   5.045

(several other examples run for ~5s on my laptop, so they may
appear in the list of elapsed time > 5s depending on cpu load)

## REQUESTED CHANGES MADE
* we changed DESCRIPTION file (packages and software written in
single quotes and names and years were added to references)