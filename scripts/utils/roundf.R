# round AND dormat to a set number of digits

roundf <- function(x, n = 0) {
  # x : real, value to round
  # n : integer, number of decimal places to return
format(round(x, digits = n), nsmall = n)
}
