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
    '        <ct:Assign>\n', x, '\n</ct:Assign>\n',
    '      </ct:Variable>\n',
    '    </mdef:StructuralModel>\n',
    '  </mdef:ModelDefinition>\n',
    '</PharmML>\n')
}

# Wrap ModelDefinition *children* in the smallest valid document.
.pharmmlWrapMdef <- function(x) {
  .pharmmlWrapMdefRaw(paste0("  <mdef:ModelDefinition>\n", x,
                             "\n  </mdef:ModelDefinition>"))
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
    '  <IndependentVariable symbId="t"/>\n', x, '\n</PharmML>\n')
}
