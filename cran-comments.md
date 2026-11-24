## Resubmission

This replaces the 0.1.0 upload made shortly before (not yet reviewed). Changes:

* Added executable examples to every exported function.
* Added URL, BugReports and the copyright-holder role.
* Reference tests now use tolerances suited to other platforms.

## Submission

This is a new submission. The package implements the method of
Xu, Jia, Chen and Chen, "Finite-Smoothing Efficiency and Inference in Single-Index Models"
(arXiv:2611.22384).

## Test environments

* Windows 11, R 4.6.2 (local)

## R CMD check results

`R CMD check --as-cran` of the first upload: 0 errors | 0 warnings | 2 notes.
For this resubmission the examples, tests, vignette and Rd checks were rerun locally and pass.

* This is a new submission.
* README.md and NEWS.md could not be checked locally because pandoc is not installed
  on the test machine.
