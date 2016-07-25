## Test environments
* local OS X install, R 3.3.1
* Ubuntu 12.04 LTS, R 3.3.1 (on travis-ci)
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs.
Building using win-builder yields the following NOTE:

	Possibly mis-spelled words in DESCRIPTION:
	  Biomedical (3:49)
	  Datasets (3:60)
	  biomedical (8:5)
	  datasets (8:16, 9:5, 9:59)

These words are not mis-spelled.

## Downstream dependencies
There are not downstream dependencies at the moment.
