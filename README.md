# latticetm
by Oliver Adams (oliver.adams@gmail.com)

This is the implementation of our [EMNLP 2016] (emnlp2016.net) paper titled:
[Learning a Lexicon and Translation Model from Phoneme Lattices] (http://people.eng.unimelb.edu.au/tcohn/papers/adams16emnlp.pdf)

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

This program builds on the codebase of [latticelm]
(https://github.com/neubig/latticelm-v2) by Graham Neubig in order to
perform translation modeling.

Install
-------

First, in terms of standard libraries, you must have autotools, libtool, and Boost. If
you are on Ubuntu/Debian linux, you can install them below:

    $ sudo apt-get install autotools libtool libboost-all

You must install OpenFST separately.

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
	$ --prior pmp --starters 0.00001 --lambda 0.75 --seed 4 \
	$ --outfile data/out/transcription

The program will run and a probabilistic transcription will be output to
data/out/transcription. Each line will be a sequence of phonemes with a space
between each phoneme. The output is probably correct, because I'm cheeky and am
using a magic seed for this example.

The lattice file (in this case data/german.lat) has `n` lattices, where
`n=train_len`. Each line specifies an arc in the form `<from> <to> <in> <out>
<prob>`. Probabilities are negative log probabilities. Blank lines delimit the
lattices. The translation file (data/english.txt) is a list of translations
corresponding to lattices.

Unfortunately we don't have license to share the BTEC data used in results
reported in the paper. In the coming months I will be applying this method to
other languages and other data sets, so I look forward to including recipes for
reproduceability. In doing so, I will endeavour to make the code more
understandable too, so bear with me.
