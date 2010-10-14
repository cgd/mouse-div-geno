#!/bin/bash

# exit on error and don't allow the use of unset variables
set -o errexit
set -o nounset

R CMD Sweave MouseDivGeno.Rnw
texi2pdf MouseDivGeno.tex
rm -f _region_.log _region_.tex MouseDivGeno.log MouseDivGeno.tex MouseDivGeno.aux MouseDivGeno.dvi MouseDivGeno.log
