function visualize(fst)
    run(`fstdraw --isymbols=symbols.txt --osymbols=symbols.txt $fst.fst $fst.dot`)
    run(pipeline(`dot -Tpdf $fst.dot`, stdout="$fst.pdf"))
    run(`evince $fst.pdf`)
end
