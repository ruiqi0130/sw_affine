# Smith-waterman based Local DNA Alignment with Affine Gap Penalties

## Installation

Install using the command below.

```R
library(devtools)
install_github("ruiqi0130/sw_affine")
```

## How to use

```R
Rscript --vanilla hw1.R <input file> <score file>
```

### Necessary input files

1. inputFile: 2 lines txt file corresponding to two sequences to be aligned.
2. scoreFile: matrix txt file, a amino-acid scoring matrix

Two additional optional inputs may be included for the opening gap and gap extension terms. Defaults for these are: `openGap=-2`, `extGap=-1`.

### Outputs

1. Two input sequences.
2. Score Matrix used for alignment
3. The best local alignment score
4. The alignment results
