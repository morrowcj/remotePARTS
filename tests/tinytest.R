# run tiny tests

if ( requireNamespace("tinytest", quietly=TRUE) ){
  tinytest::test_package("remotePARTS")
}
