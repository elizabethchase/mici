resample <- function(x, ...) x[sample.int(length(x), ...)]

`%notin%` <- Negate(`%in%`)