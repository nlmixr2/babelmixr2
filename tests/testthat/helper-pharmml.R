# Shared helpers for the PharmML tests.  testthat sources helper-*.R before the
# test files, so these are visible to all of them.

# Wrap a math subtree in the smallest document the schema will accept, so the
# math walker's output can be validated in isolation.
.pharmmlWrapMath <- function(x) {
  paste0(
    '<?xml version="1.0" encoding="UTF-8"?>\n',
    '<PharmML xmlns="http://www.pharmml.org/pharmml/0.9/PharmML"\n',
    '    xmlns:ct="http://www.pharmml.org/pharmml/0.9/CommonTypes"\n',
    '    xmlns:math="http://www.pharmml.org/pharmml/0.9/Maths"\n',
    '    xmlns:mdef="http://www.pharmml.org/pharmml/0.9/ModelDefinition"\n',
    '    writtenVersion="0.9" id="i1">\n',
    '  <ct:Name>math walker fixture</ct:Name>\n',
    '  <IndependentVariable symbId="t"/>\n',
    '  <mdef:ModelDefinition>\n',
    '    <mdef:StructuralModel blkId="sm1">\n',
    '      <ct:Variable symbId="y" symbolType="real">\n',
    '        <ct:Assign>\n',
    x,
    '\n</ct:Assign>\n',
    '      </ct:Variable>\n',
    '    </mdef:StructuralModel>\n',
    '  </mdef:ModelDefinition>\n',
    '</PharmML>\n'
  )
}

# Wrap ModelDefinition *children* in the smallest valid document.
.pharmmlWrapMdef <- function(x) {
  .pharmmlWrapMdefRaw(paste0(
    "  <mdef:ModelDefinition>\n",
    x,
    "\n  </mdef:ModelDefinition>"
  ))
}

# Wrap an already-complete <mdef:ModelDefinition> element.
.pharmmlWrapMdefRaw <- function(x) {
  paste0(
    '<?xml version="1.0" encoding="UTF-8"?>\n',
    '<PharmML xmlns="http://www.pharmml.org/pharmml/0.9/PharmML"\n',
    '    xmlns:ct="http://www.pharmml.org/pharmml/0.9/CommonTypes"\n',
    '    xmlns:math="http://www.pharmml.org/pharmml/0.9/Maths"\n',
    '    xmlns:mdef="http://www.pharmml.org/pharmml/0.9/ModelDefinition"\n',
    '    xmlns:po="http://www.pharmml.org/probonto/ProbOnto"\n',
    '    writtenVersion="0.9" id="i1">\n',
    '  <ct:Name>fixture</ct:Name>\n',
    '  <IndependentVariable symbId="t"/>\n',
    x,
    '\n</PharmML>\n'
  )
}

# Every `SymbRef` that names a block must name a symbol that block declares.
# Returns the dangling references as "blk/symb", so a test can expect none.
.pharmmlDanglingRefs <- function(doc) {
  .x <- xml2::read_xml(as.character(doc))
  .refs <- xml2::xml_find_all(.x, "//*[local-name()='SymbRef'][@blkIdRef]")
  .ret <- vapply(
    .refs,
    function(.r) {
      .blk <- xml2::xml_attr(.r, "blkIdRef")
      .sym <- xml2::xml_attr(.r, "symbIdRef")
      .found <- xml2::xml_find_all(
        .x,
        sprintf("//*[@blkId='%s']//*[@symbId='%s']", .blk, .sym)
      )
      if (length(.found) > 0L) {
        return(NA_character_)
      }
      paste0(.blk, "/", .sym)
    },
    character(1)
  )
  unique(.ret[!is.na(.ret)])
}
