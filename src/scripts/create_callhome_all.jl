#!/home/oadams/code/julia/julia

corpus_dir="data/fisher-callhome/corpus"
test_num=1818

function create_callhome_all()
    run(pipeline(`cat $corpus_dir/plf/callhome_train.es $corpus_dir/plf/callhome_devtest.es $corpus_dir/plf/callhome_evltest.es`, stdout="$corpus_dir/plf/callhome_all.es"))
    run(pipeline(`cat $corpus_dir/ldc/callhome_train.es $corpus_dir/ldc/callhome_devtest.es $corpus_dir/ldc/callhome_evltest.es`, stdout="$corpus_dir/ldc/callhome_all.es"))
    run(pipeline(`cat $corpus_dir/ldc/callhome_train.en $corpus_dir/ldc/callhome_devtest.en $corpus_dir/plf/callhome_evltest.en`, stdout="$corpus_dir/ldc/callhome_all.en"))
    run(`./plf2txt.py $corpus_dir callhome_all`)
end

function create_callhome_subsets()
    for i = 0:10
        n = test_num+i*test_num
        run(pipeline(`tail -n $n $corpus_dir/ldc/callhome_all.es.nopunc.lower`, stdout="$corpus_dir/ldc/callhome_all.es.nopunc.lower.n$n"))
        run(pipeline(`tail -n $n $corpus_dir/ldc/callhome_all.en.nopunc.lower`, stdout="$corpus_dir/ldc/callhome_all.en.nopunc.lower.n$n"))
        run(pipeline(`tail -n $n $corpus_dir/plf/callhome_all.es.lat`, stdout="$corpus_dir/ldc/callhome_all.es.lat.n$n"))
    end
end

function train_giza_subsets()
    for i = 0:10
        n = test_num+i*test_num
        run(`./run_giza callhome_all $n`)
    end
end

#train_giza_subsets()

function decode_evaluate_giza_subsets()
    for i = 0:10
        n = test_num+i*test_num
        tm_file = "data/out/giza-pp/callhome_all/n$n/train/model/lex.e2f"
        cmd=pipeline(`./decode_giza callhome_evltest $tm_file`,
                stdout="data/out/giza-pp/callhome_all/n$n/callhome_evltest.es")
        println(cmd)
        run(cmd)
        cmd=`python ../nlp/per.py --ref data/fisher-callhome/corpus/ldc/callhome_evltest.es.nopunc.lower --hypo data/out/giza-pp/callhome_all/n$n/callhome_evltest.es`
        println(cmd)
        run(cmd)
    end
end

decode_evaluate_giza_subsets()
