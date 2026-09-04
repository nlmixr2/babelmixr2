test_that(".pmlNode emits empty, attributed and nested nodes", {
  expect_equal(.pmlNode("ct:Real"), "<ct:Real/>")

  expect_equal(.pmlNode("ct:SymbRef", attrs = c(symbIdRef = "V")),
               '<ct:SymbRef symbIdRef="V"/>')

  expect_equal(.pmlNode("ct:SymbRef", attrs = c(blkIdRef = "pm1", symbIdRef = "V")),
               '<ct:SymbRef blkIdRef="pm1" symbIdRef="V"/>')

  expect_equal(.pmlNode("ct:Assign", children = "<ct:Real>1</ct:Real>"),
               "<ct:Assign>\n<ct:Real>1</ct:Real>\n</ct:Assign>")
})

test_that(".pmlNode indents", {
  expect_equal(.pmlNode("ct:Real", indent = 2L), "        <ct:Real/>")
})

test_that(".pmlNode escapes attribute values", {
  expect_equal(.pmlNode("x", attrs = c(a = 'q"&<>')),
               '<x a="q&quot;&amp;&lt;&gt;"/>')
})

test_that(".pmlText emits a text node", {
  expect_equal(.pmlText("ct:Real", 1.5), "<ct:Real>1.5</ct:Real>")
  expect_equal(.pmlText("ct:Int", 3L), "<ct:Int>3</ct:Int>")
})

test_that(".rxToPharmml handles scalars", {
  expect_equal(.rxToPharmml(quote(1L)), "<ct:Int>1</ct:Int>")
  expect_equal(.rxToPharmml(quote(1.5)), "<ct:Real>1.5</ct:Real>")
  expect_equal(.rxToPharmml(2), "<ct:Real>2</ct:Real>")
})

test_that(".rxToPharmml handles symbols", {
  expect_equal(.rxToPharmml(quote(V)), '<ct:SymbRef symbIdRef="V"/>')
})

test_that(".rxToPharmml maps time to the independent variable", {
  expect_equal(.rxToPharmml(quote(time)), '<ct:SymbRef symbIdRef="t"/>')
})

test_that(".rxToPharmml inlines numeric constants", {
  expect_equal(.rxToPharmml(quote(pi)), "<ct:Real>3.141592653589793</ct:Real>")
  expect_equal(.rxToPharmml(quote(M_LN2)), "<ct:Real>0.6931471805599453</ct:Real>")
})

test_that(".rxToPharmml rejects constructs PharmML cannot express", {
  expect_error(.rxToPharmml(quote(linCmt())), "linCmt")
  expect_error(.rxToPharmml(quote(NA)), "NA")
  expect_error(.rxToPharmml(quote(digamma(x))), "digamma")
})

test_that(".rxToPharmml handles binary operators", {
  expect_equal(
    .rxToPharmml(quote(a + b)),
    paste0("<math:Binop op=\"plus\">\n",
           "<ct:SymbRef symbIdRef=\"a\"/>\n",
           "<ct:SymbRef symbIdRef=\"b\"/>\n",
           "</math:Binop>"))

  expect_match(.rxToPharmml(quote(a * b)), 'op="times"')
  expect_match(.rxToPharmml(quote(a / b)), 'op="divide"')
  expect_match(.rxToPharmml(quote(a - b)), 'op="minus"')
  expect_match(.rxToPharmml(quote(a ^ b)), 'op="power"')
})

test_that(".rxToPharmml handles unary minus", {
  expect_equal(
    .rxToPharmml(quote(-a)),
    paste0("<math:Uniop op=\"minus\">\n",
           "<ct:SymbRef symbIdRef=\"a\"/>\n",
           "</math:Uniop>"))
})

test_that(".rxToPharmml maps rxode2 functions onto native Uniop values", {
  expect_match(.rxToPharmml(quote(exp(a))),   'Uniop op="exp"')
  expect_match(.rxToPharmml(quote(log(a))),   'Uniop op="log"')
  expect_match(.rxToPharmml(quote(sqrt(a))),  'Uniop op="sqrt"')
  expect_match(.rxToPharmml(quote(log10(a))), 'Uniop op="log10"')
  # PharmML has these natively -- unlike Monolix, no rewrite needed
  expect_match(.rxToPharmml(quote(logit(a))),  'Uniop op="logit"')
  expect_match(.rxToPharmml(quote(expit(a))),  'Uniop op="logistic"')
  expect_match(.rxToPharmml(quote(probit(a))), 'Uniop op="probit"')
  expect_match(.rxToPharmml(quote(pnorm(a))),  'Uniop op="normcdf"')
  expect_match(.rxToPharmml(quote(lgamma(a))), 'Uniop op="gammaln"')
})

test_that(".rxToPharmml rewrites functions with no direct Uniop", {
  .x <- .rxToPharmml(quote(log1p(a)))
  expect_match(.x, 'Uniop op="log"')
  expect_match(.x, 'Binop op="plus"')

  .x <- .rxToPharmml(quote(expm1(a)))
  expect_match(.x, 'Uniop op="exp"')
  expect_match(.x, 'Binop op="minus"')
})

test_that(".rxToPharmml handles two-argument functions", {
  expect_match(.rxToPharmml(quote(atan2(a, b))), 'Binop op="atan2"')
  expect_match(.rxToPharmml(quote(max(a, b))),   'Binop op="max"')
  expect_match(.rxToPharmml(quote(min(a, b))),   'Binop op="min"')
})

test_that(".rxToPharmml handles parenthesised expressions transparently", {
  expect_equal(.rxToPharmml(quote((a))), '<ct:SymbRef symbIdRef="a"/>')
  expect_match(.rxToPharmml(quote((a + b) * c)), 'op="times"')
})

test_that(".rxToPharmml nests correctly and preserves precedence", {
  .x <- .rxToPharmml(quote(ka * depot - cl / v * central))
  expect_match(.x, "^<math:Binop op=\"minus\">")
  expect_equal(lengths(regmatches(.x, gregexpr("math:Binop", .x))), 8L)
})

test_that(".rxToPharmml maps logical operators", {
  expect_match(.rxToPharmml(quote(a < b)),  'LogicBinop op="lt"')
  expect_match(.rxToPharmml(quote(a <= b)), 'LogicBinop op="leq"')
  expect_match(.rxToPharmml(quote(a > b)),  'LogicBinop op="gt"')
  expect_match(.rxToPharmml(quote(a >= b)), 'LogicBinop op="geq"')
  expect_match(.rxToPharmml(quote(a == b)), 'LogicBinop op="eq"')
  expect_match(.rxToPharmml(quote(a != b)), 'LogicBinop op="neq"')
  expect_match(.rxToPharmml(quote(a & b)),  'LogicBinop op="and"')
  expect_match(.rxToPharmml(quote(a | b)),  'LogicBinop op="or"')
  expect_match(.rxToPharmml(quote(!a)),     'LogicUniop op="not"')
})

test_that(".rxToPharmml maps ifelse() to a two-piece Piecewise", {
  .x <- .rxToPharmml(quote(ifelse(wt > 70, a, b)))
  expect_match(.x, "^<math:Piecewise>")
  expect_equal(lengths(regmatches(.x, gregexpr("<math:Piece>", .x))), 2L)
  expect_match(.x, "<math:Otherwise/>")
  expect_match(.x, 'LogicBinop op="gt"')
})

test_that(".rxToPharmml maps a bare if/else to Piecewise", {
  .x <- .rxToPharmml(quote(if (wt > 70) a else b))
  expect_equal(lengths(regmatches(.x, gregexpr("<math:Piece>", .x))), 2L)
  expect_match(.x, "<math:Otherwise/>")
})

test_that(".rxToPharmml maps an if with no else to a one-piece Piecewise", {
  .x <- .rxToPharmml(quote(if (wt > 70) a))
  expect_equal(lengths(regmatches(.x, gregexpr("<math:Piece>", .x))), 1L)
  expect_false(grepl("Otherwise", .x, fixed = TRUE))
})

test_that(".rxToPharmml handles nested ifelse()", {
  .x <- .rxToPharmml(quote(ifelse(a > 1, x, ifelse(b > 2, y, z))))
  expect_equal(lengths(regmatches(.x, gregexpr("<math:Piecewise>", .x))), 2L)
})

# Wrap a math subtree in the smallest document the schema will accept, so the
# walker's output can be validated in isolation.
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

test_that("the math walker emits schema-valid PharmML", {
  .exprs <- list(
    quote(a + b),
    quote(-a),
    quote(ka * depot - cl / v * central),
    quote(exp(tcl + eta.cl)),
    quote(logit(a)),
    quote(expit(a)),
    quote(log1p(a)),
    quote(atan2(a, b)),
    quote(ifelse(wt > 70, a, b)),
    quote(if (wt > 70) a else b),
    quote(ifelse(a > 1, x, ifelse(b > 2, y, z)))
  )
  for (.e in .exprs) {
    .doc <- .pharmmlWrapMath(.rxToPharmml(.e))
    expect_true(pharmmlValidate(.doc),
                info = paste("failed for:", deparse1(.e)))
  }
})
