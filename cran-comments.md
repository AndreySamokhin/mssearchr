## Resubmission

This resubmission addresses several recommendations from the CRAN team.

* It was recommended to include a web reference for the API in the DESCRIPTION
file, but there is no official webpage dedicated to the API. The API details are
solely available in the software manual. Therefore, a web reference to 'MS
Search' software (NIST) has been added to the DESCRIPTION file.
* It was recommended to include references describing the methods used in the
package. Currently, there are no references available for the package as an
article detailing these methods is currently being written.
* It was recommended to replace '\dontrun' with '\donttest' where possible.
Currently, '\dontrun' is used only for the 'LibrarySearchUsingNistApi()'
function because the respective example cannot be executed without MS Search
(NIST) software installed.



## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.
