
monolixGetData <- function(ui, data) {
  dat <- rxode2::etTrans(data, ui$simulationModel, addCmt=TRUE, allTimeVar=TRUE)
}
