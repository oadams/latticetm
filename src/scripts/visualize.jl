function visualize_tm()
    fst="tm"
    run(`fstdraw --isymbols=f_symbols.txt --osymbols=e_symbols.txt $fst.fst $fst.dot`)
    run(pipeline(`dot -Tpdf $fst.dot`, stdout="$fst.pdf"))
    run(`evince $fst.pdf`)
end

function visualize_latlex()
    fst="latlex"
    run(`fstdraw --isymbols=f_symbols.txt --osymbols=f_symbols.txt $fst.fst $fst.dot`)
    run(pipeline(`dot -Tpdf $fst.dot`, stdout="$fst.pdf"))
    run(`evince $fst.pdf`)
end

function visualize_lexicon()
    fst="lexicon"
    run(`fstdraw --isymbols=f_symbols.txt --osymbols=f_symbols.txt $fst.fst $fst.dot`)
    run(pipeline(`dot -Tpdf $fst.dot`, stdout="$fst.pdf"))
    run(`evince $fst.pdf`)
end

function visualize_composed()
    fst="composed"
    run(`fstdraw --isymbols=f_symbols.txt --osymbols=e_symbols.txt $fst.fst $fst.dot`)
    run(pipeline(`dot -Tpdf $fst.dot`, stdout="$fst.pdf"))
    run(`evince $fst.pdf`)
end

function visualize_sample()
    fst="sample"
    run(`fstdraw --isymbols=f_symbols.txt --osymbols=e_symbols.txt $fst.fst $fst.dot`)
    run(pipeline(`dot -Tpdf $fst.dot`, stdout="$fst.pdf"))
    run(`evince $fst.pdf`)
end
