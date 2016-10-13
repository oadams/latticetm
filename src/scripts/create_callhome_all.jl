#!/home/oadams/code/julia/julia

corpus_dir="data/fisher-callhome/corpus"
callhome_test_num=1818
fisher_test_num=3629

function create_callhome_all()
    run(pipeline(`cat $corpus_dir/plf/callhome_train.es $corpus_dir/plf/callhome_devtest.es $corpus_dir/plf/callhome_evltest.es`, stdout="$corpus_dir/plf/callhome_all.es"))
    run(pipeline(`cat $corpus_dir/ldc/callhome_train.es $corpus_dir/ldc/callhome_devtest.es $corpus_dir/ldc/callhome_evltest.es`, stdout="$corpus_dir/ldc/callhome_all.es"))
    run(pipeline(`cat $corpus_dir/ldc/callhome_train.en $corpus_dir/ldc/callhome_devtest.en $corpus_dir/plf/callhome_evltest.en`, stdout="$corpus_dir/ldc/callhome_all.en"))
    run(`./plf2txt.py $corpus_dir callhome_all`)
end

function create_fisher_all()
    run(pipeline(`cat $corpus_dir/plf/fisher_train.es $corpus_dir/plf/fisher_dev.es $corpus_dir/plf/fisher_dev2.es $corpus_dir/plf/fisher_test.es`, stdout="$corpus_dir/plf/fisher_all.es"))
    run(pipeline(`cat $corpus_dir/ldc/fisher_train.es $corpus_dir/ldc/fisher_dev.es $corpus_dir/ldc/fisher_dev2.es $corpus_dir/ldc/fisher_test.es`, stdout="$corpus_dir/ldc/fisher_all.es"))
    run(pipeline(`cat $corpus_dir/ldc/fisher_train.en $corpus_dir/ldc/fisher_dev.en $corpus_dir/ldc/fisher_dev2.en $corpus_dir/ldc/fisher_test.en`, stdout="$corpus_dir/ldc/fisher_all.en"))
    run(`src/scripts/plf2txt.py $corpus_dir fisher_all`)
end

function create_callhome_subsets()
    for i = 0:10
        n = callhome_test_num+i*callhome_test_num
        run(pipeline(`tail -n $n $corpus_dir/ldc/callhome_all.es.nopunc.lower`, stdout="$corpus_dir/ldc/callhome_all.es.nopunc.lower.n$n"))
        run(pipeline(`tail -n $n $corpus_dir/ldc/callhome_all.en.nopunc.lower`, stdout="$corpus_dir/ldc/callhome_all.en.nopunc.lower.n$n"))
        run(pipeline(`tail -n $n $corpus_dir/plf/callhome_all.es.lat`, stdout="$corpus_dir/ldc/callhome_all.es.lat.n$n"))
    end
end

function create_fisher_subsets()
    for i = 0:10
        n = fisher_test_num+i*fisher_test_num
        run(pipeline(`tail -n $n $corpus_dir/ldc/fisher_all.es.nopunc.lower`, stdout="$corpus_dir/ldc/fisher_all.es.nopunc.lower.n$n"))
        run(pipeline(`tail -n $n $corpus_dir/ldc/fisher_all.en.nopunc.lower`, stdout="$corpus_dir/ldc/fisher_all.en.nopunc.lower.n$n"))
        run(pipeline(`tail -n $n $corpus_dir/plf/fisher_all.es.lat`, stdout="$corpus_dir/ldc/fisher_all.es.lat.n$n"))
    end
end

function create_fisher_50k_subsets()
    for i = 0:10
        n = fisher_test_num+i*fisher_test_num
        run(pipeline(`tail -n $n $corpus_dir/ldc/fisher_50k.es.nopunc.lower`, stdout="$corpus_dir/ldc/fisher_50k.es.nopunc.lower.n$n"))
        run(pipeline(`tail -n $n $corpus_dir/ldc/fisher_50k.en.nopunc.lower`, stdout="$corpus_dir/ldc/fisher_50k.en.nopunc.lower.n$n"))
        run(pipeline(`tail -n $n $corpus_dir/plf/fisher_50k.es.lat`, stdout="$corpus_dir/ldc/fisher_50k.es.lat.n$n"))
    end
end

HOME="/home/oadams"
giza_dir="$HOME/tools/giza-pp"
moses_dir="$HOME/tools/mosesdecoder"
num_cores=3

function train_giza(dataset, num_sents)
    src="data/out/plain_best_paths/$dataset.es"
    tgt="data/fisher-callhome/corpus/ldc/$dataset.en.nopunc.lower"
    OUT_DIR="$HOME/code/latticelm-v2/data/out/giza-pp/$dataset/n$num_sents"
    run(`mkdir -p $OUT_DIR/train`)
    run(`touch $OUT_DIR/fakelm`)
    run(pipeline(`tail -n $num_sents $src`, stdout="$OUT_DIR/$dataset.es"))
    run(pipeline(`tail -n $num_sents $tgt`, stdout="$OUT_DIR/$dataset.en"))
    run(`$moses_dir/scripts/training/train-model.perl \
        -root-dir $OUT_DIR/train \
        -cores $num_cores \
        -corpus $OUT_DIR/$dataset -f es -e en \
        -alignment grow-diag-final-and -reordering msd-bidirectional-fe \
        -lm 0:5:$OUT_DIR/fakelm \
        -external-bin-dir $giza_dir/bin &> $OUT_DIR/training.out`)
end

function train_giza_subsets(dataset, test_num)
    for i = 0:10
        n = test_num+i*test_num
        train_giza(dataset, n)
    end
end

function decode_evaluate_giza_subsets(trainset, evlset, test_num)
    #sents_per_hour=20875/(20)
    sents_per_hour=150400/(163)
    for i = 0:10
        n = test_num+i*test_num
        tm_file = "data/out/giza-pp/$trainset/n$n/train/model/lex.e2f"
        cmd=pipeline(`./decode_giza $evlset $tm_file`,
                stdout="data/out/giza-pp/$trainset/n$n/$evlset.es")
        #println(cmd)
        run(cmd)
        cmd=`python ../nlp/per.py --ref data/fisher-callhome/corpus/ldc/$evlset.es.nopunc.lower --hypo data/out/giza-pp/$trainset/n$n/$evlset.es`
        #println(cmd)
        wer=chomp(readall(cmd))
        hours=n/sents_per_hour
        println("(",hours,",",wer,")")
    end
end

function latticetm_callhome_scaleup(train_set, test_set, test_num)
    sents_per_hour=20875/(20)
    epochs=30
    alpha=1
    weight=1
    for i = 0:10
        train_num = test_num+i*test_num
        run(`./run $train_set $epochs $alpha $weight $train_num $test_num`)
        cmd=`./wer $train_set $epochs $alpha $weight $train_num $test_set $test_num`
        wer=chomp(readall(cmd))
        hours=train_num/sents_per_hour
        println("(",hours,",",wer,")")
    end
end

function latticetm_fisher_scaleup(train_set, test_set, test_num)
    sents_per_hour=150400/(163)
    epochs=30
    alpha=1
    weight=1
    for i = 0:10
        train_num = test_num+i*test_num
        run(`./run $train_set $epochs $alpha $weight $train_num $test_num`)
        cmd=`./wer $train_set $epochs $alpha $weight $train_num $test_set $test_num`
        wer=chomp(readall(cmd))
        hours=train_num/sents_per_hour
        println("(",hours,",",wer,")")
    end
end

#train_giza_subsets("callhome_all")
#decode_evaluate_giza_subsets("callhome_all", "callhome_evltest")
#latticetm_scaleup("callhome_all", "callhome_evltest")

#create_fisher_all()
#create_fisher_50k_subsets()
train_giza_subsets("fisher_50k", fisher_test_num)
#decode_evaluate_giza_subsets("fisher_50k", "fisher_test", fisher_test_num)
#latticetm_fisher_scaleup("fisher_50k", "fisher_test", fisher_test_num)
