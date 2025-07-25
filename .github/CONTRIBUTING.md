# Contributing to remotePARTS

Here, we outline how to suggest changes to remotePARTS. 

## Issues

For all suggestions, including bug reports, feature requests, documentation, 
and typo or language corrections, please submit an issue on Github. 

## Code contribution

Substantial changes can be suggested by submitting a pull request into the
`develop` branch, provided that there is an associated (linked) issue. If 
there is not an existing issue concerning the changes provided in the pull 
request, you should create the issue first and link it in the pull request. 

Ideally, pull requests should be small in scope (i.e., dedicated to a single
issue/feature/functionality).

## Tests

Before submitting code (via pull request), please check that 1) your branch
is up-to-date with the destination branch, and 2) automated tests were 
successful. Use `devtools::test()` to run unit tests and `devtools::check()` 
to check the full package. The main branches will also run these checks
automatically, via GitHub actions, when pull requests are submitted.
