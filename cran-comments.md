## Resubmission

This fixes the problems found by the incoming pre-tests of the previous 0.1.0 upload:

* A stray control character in index_predict.Rd broke the PDF and HTML manuals; it is
  replaced by the intended `\beta`. The PDF manual now builds without errors.

The flagged words "Jia" and "Xu" in DESCRIPTION are author names.

## Test environments

* Windows 11, R 4.6.2 (local), with MiKTeX for the PDF manual

## R CMD check results

0 errors | 0 warnings | 3 notes

* This is a new submission.
* README.md and NEWS.md could not be checked locally (pandoc is not installed locally).
* HTML validation was skipped locally (HTML Tidy is not installed locally).
