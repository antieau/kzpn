# KZPN
`SAGE` code and example worksheet for computing the p-adic syntomic cohomology and p-adic algebraic
K-groups of quotients of rings of integers in totally ramified p-adic fields.

This code has been verified to work in `SAGE 10.3`.



## Instructions

### Installation

Ensure that `SAGE` is
[installed](https://doc.sagemath.org/html/en/installation/index.html).
Clone this `git` repository by running the following command:

```
git clone git@github.com:antieau/kzpn.git
```

### Testing the installation

In the `\kzpn\` folder, run

```
sage test_main.py
```

Several tests are run and they should each result in a `PASS`.


### Using

There are three examples in the main directory to illustrate how to use the
code. They can be run by copying them into an `iPython` notebook (running with
the `SAGE` kernel) or by running, for example,

```
sage example_2101_n2_i2.py
```

from the command-line.



## Bibliography

Antieau, Krause, Nikolaus, _On the K-theory of $$\mathbb{Z}/p^n$$ -- announcement_, \[[arXiv:2204](https://arxiv.org/abs/2204.03420)\].

Antieau, Krause, Nikolaus, _On the K-theory of $$\mathbb{Z}/p^n$$_, \[[arXiv:2405]()\].

Antieau, Krause, Nikolaus, _Prismatic cohomology relative to $$\delta$$-rings_, \[[arXiv:2310](https://arxiv.org/abs/2310.12770)\].

SageMath, the Sage Mathematics Software System (Version 10.3), The Sage Developers, 2024, https://www.sagemath.org.
