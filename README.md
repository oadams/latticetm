# latticetm
by Oliver Adams (oliver.adams@gmail.com)

This is the implementation of our [EMNLP 2016] (emnlp2016.net) paper titled: 
*Learning a Lexicon and Translation Model from Phoneme Lattices*.

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
additionally perform translation modeling as well as segmentation of phonemes.

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
