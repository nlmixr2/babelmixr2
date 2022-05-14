.getUiFunFromIniAndModel <- function(ui, ini, model) {
  .ls <- ls(ui$meta, all.names=TRUE)
  .ret <- vector("list", length(.ls) + 3)
  .ret[[1]] <- quote(`{`)
  for (.i in seq_along(.ls)) {
    .ret[[.i + 1]] <- eval(parse(text=paste("quote(", .ls[.i], "<-", deparse1(ui$meta[[.ls[.i]]]), ")")))
  }
  .len <- length(.ls)
  .ret[[.len + 2]] <- ini
  .ret[[.len + 3]] <- model
  .retf <- function(){}
  body(.retf) <- as.call(.ret)
  .retf
}
