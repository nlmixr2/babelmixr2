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
