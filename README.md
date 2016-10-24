# latticetm
by Oliver Adams (oliver.adams@gmail.com), building on the codebase of [latticelm]
(https://github.com/neubig/latticelm-v2) by Graham Neubig.

This is the implementation of our [EMNLP 2016] (emnlp2016.net) paper,
[Learning a Lexicon and Translation Model from Phoneme Lattices] (http://people.eng.unimelb.edu.au/tcohn/papers/adams16emnlp.pdf), which won the best short paper award.

If you use this code, please cite the paper

```
@inproceedings{adams16emnlp,
    title = {Learning a Lexicon and Translation Model from Phoneme Lattices},
    author = {Oliver Adams and Graham Neubig and Trevor Cohn and Steven Bird and Quoc Truong Do and Satoshi Nakamura},
    booktitle = {Conference on Empirical Methods in Natural Language Processing (EMNLP)},
    address = {Austin, Texas, USA},
    month = {November},
    year = {2016}
}
```

This program builds on the codebase of [latticelm] (https://github.com/neubig/latticelm-v2) in order to
perform translation modeling.

Install
-------

First, in terms of standard libraries, you must have autotools, libtool, and Boost. If
you are on Ubuntu/Debian linux, you can install them below:

    $ sudo apt-get install autotools libtool libboost-all-dev

You must install [OpenFST] (http://www.openfst.org/) separately.

Once these two packages are installed, run the following commands, specifying the
correct path for openfst (likely /usr/local/ on Debian-based systems).

    $ autoreconf -i
    $ ./configure --with-openfst=/path/to/openfst
    $ make

Usage
-----

A toy dataset is available in data/. It is the same example as in the paper,
and illustrates the formatting of input files to the program. To run this
example:

	$ ./src/latticelm/latticelm \
	$ --train_file data/german.lat --trans_file data/english.txt \
	$ --file_format openfst --model_type lextm \
	$ --epochs 11 --concentration 1 --lattice_weight 1 \
	$ --train_len 3 --test_len 3 \
	$ --prior pmp --starters 0.00001 --gamma 0.75 --seed 4 \
	$ --outfile data/out/transcription

The program will run and a probabilistic transcription will be output to
data/out/transcription. Each line will be a sequence of phonemes with a space
between each phoneme (segmentation isn't currently shown in the transcript).
The output is probably correct, because I'm cheeky and am using a magic seed
for this example.

The lattice file (in this case data/german.lat) has `n` lattices, where
`n=train_len`. Each line specifies an arc in the form `<from> <to> <in> <out>
<prob>` (this is referred to as the `openfst` format). Probabilities are negative log probabilities. Blank lines delimit the
lattices. The translation file (data/english.txt) is a list of translations
corresponding to lattices.

Other arguments include:
- `--epochs` is the number of iterations of the corpus for sampling.
- `--concentration` is the Dirichlet process concentration parameter (ie. alpha in the paper, giving the strength of the base distribution).
- `--lattice_weight` tweaks how much importance is given to lattice weights. Don't worry about it, just leave it at 1.
- `--test_len` specifies how many lines you actually want transcribed. This is useful for keeping a uniform test set that is a subset of some larger corpora. For example, we can have 1,000 test set lines, but train on larger supersets.
- `--prior` is the spelling model prior. options are `geom` for geometric, `poisson` for poisson, and `pmp` for what is called *shifted* in the paper.
- `--starters` is a hyperparam relevant only to `pmp` (*shifted*) prior. It is a list of k floats that specify the base probability of a word of length 1..k respectively. In the example above, there is only one probability specified to decrease the chance of words of length one. More probabilities can be entered such as `--starters 0.00001 0.1` if you want the probability of a word of length 2 to be 0.1. The remaining probability mass is distributed geometrically. `"pmp"` is an acronym for poor-man's Poisson. It is actually faster though and frequently outperforms the Poisson spelling model.
- `--lambda` is the `poisson` hyperparam. Not included in the above example, because it uses `pmp`.
- `--gamma` is the hyperparam that describes decay for the geometric (`geom`) distribution and the *shifted* distribution (`pmp`).
- `--seed` is the seed for the random number generator so that certain results can be reproduced.
- `--outfile` is the file to put the *unsegmented* automatic transcription that harnesses the provided translation.

Unfortunately we don't have license to share the BTEC data used in results
reported in the paper. In the coming months I will be applying this method to
other languages and other data sets, so I look forward to including recipes for
reproduceability. In doing so, I will endeavour to make the code more
understandable too, so bear with me.
