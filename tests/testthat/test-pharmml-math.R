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
