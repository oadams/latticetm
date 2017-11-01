exp_dir="/home/oadams/mam/exp/573"
./src/latticelm/latticelm \
	--train_file $exp_dir/lattices/lattice_fn_list.txt \
	--trans_file $exp_dir/lattices/transl_fn_list.txt \
	--file_format openfst --model_type lextm  --epochs 11 --concentration 1 \
	--lattice_weight 1 --train_len 3 --test_len 3 --prior pmp \
	--starters 0.00001 --gamma 0.75 --seed 4 \
	--plain_best_paths best_path \
	--symbol_file $exp_dir/lattices/symbols.txt

