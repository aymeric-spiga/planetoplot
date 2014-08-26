#!/bin/bash
# This script can be used to drive the bam.sh script.
# It also acts as a description of the different parts of the bam script.
# Remember to check path to files in bam.sh.

# bam.sh boundaries #
# generic bounding condition check
./bam.sh boundaries temp # temperature of outer nest boundaries
./bam.sh boundaries uv # horizontal winds of outer nest boundaries

# bam.sh sync #
# Synchronization check using tsurf 
./bam.sh sync

# bam.sh nests #
# generic parent to child nest comparisons
./bam.sh nests power # power spectrum of horizontal winds between nests
./bam.sh nests histo # histograms of horizontal winds between nests

# bam.sh srlsn3 #
# generic analysis at the Safe Reference Landing Site (srls) for nest 3 (n3)
./bam.sh srlsn3 power # power spectrum of horizontal winds between LMD_MMM and MRAMS
./bam.sh srlsn3 histo # histograms of horizontal winds between LMD_MMM and MRAMS for each descent phase
./bam.sh srlsn3 series # time-series (1d and 2d) at landing site (winds, temp, hodograph, wind rotation)
./bam.sh srlsn3 profiles # profiles at landing site
./bam.sh srlsn3 surface # NOT TESTED YET BECAUSE OF THE LACK OF FIELDS

# bam.sh maps
# generic analysis of lat/lon maps for nest 3
./bam.sh maps winds
./bam.sh maps surface  # NOT TESTED YET BECAUSE OF THE LACK OF FIELDS

# bam.sh variability
# day-to-day variability, studied with the day before and the day after
./bam.sh variability histo # histograms of horizontal winds of the 3 days, + comp LMD/MRAMS
./bam.sh variability density # kernel density estimates of day-to-day differences => standard deviation, mean, kurtosis, skewness
