# CommutativeMonoids
This site pretends to be a site for the study of finitely generated commutative monoids, providing python-notebook for that.

At the moment the files you can find are:

* `README.md`: this file.

* `ImportFromBitbucket.py`: an small library for loading the  

* `NumericalSemigroups.ipynb`: a notebook with the class `NumericalSemigroup`.

* `NumericalSemigroups.py`: the library obtained from the above notebook.

* `integerSmithNormalFormAndApplications.ipynb`: functions for computing the integer Smith normal form that is used mainly in cancellative monoids.

* `integerSmithNormalFormAndApplications.py`: the library obtained from the above notebook.

* `minimalsNp`: Computation of the minimal elements of subsets of `\mathbb{N}^p`.

* `AffineSemigroup.ipynb`: a notebook with the class `AffineSemigroup`.

To use any of this files jush install the package `ImportFromBitbucket` from pypi with the command `pip install ImportFromBitbucket`.
Once installed the above package you can use the libraries. For instance to use the functions of `integerSmithNormalFormAndApplications`
run the commands:

* `import ImportFromBitbucket`

* `ImportFromBitbucket.loadPyFile('integerSmithNormalFormAndApplications.py')`

The contributors to this project are:

* J. I. García-García, Universidad de Cádiz, ignacio.garcia@uca.es

* D. Marín-Aragón, Universidad de Cádiz, daniel.marin@uca.es

* A. Sánchez-Roselly Navarro, Universidad de Cádiz, alfredo.sanchez@uca.es

* A. Vigneron-Tenorio, Universidad de Cádiz, alberto.vigneron@uca.es

You can visualize the above notebooks files in the following links:

* <a href='http://nbviewer.jupyter.org/urls/bitbucket.org/juan_ignacio_garcia_garcia/commutativemonoids/raw/master/NumericalSemigroups.ipynb?flush_cache=true' target='_blank'>NumericalSemigroups.ipynb</a>

* <a href='http://nbviewer.jupyter.org/urls/bitbucket.org/juan_ignacio_garcia_garcia/commutativemonoids/raw/master/integerSmithNormalFormAndApplications.ipynb?flush_cache=true' target='_blank'>integerSmithNormalFormAndApplications.ipynb</a>

* <a href='http://nbviewer.jupyter.org/urls/bitbucket.org/juan_ignacio_garcia_garcia/commutativemonoids/raw/master/minimalsNp.ipynb?flush_cache=true' target='_blank'>minimalsNp.py</a>

* <a href='http://nbviewer.jupyter.org/urls/bitbucket.org/juan_ignacio_garcia_garcia/commutativemonoids/raw/master/AffineSemigroup.ipynb?flush_cache=true' target='_blank'>AffineSemigroup.ipynb</a>

This repository is a migration of the repository https://bitbucket.org/juan_ignacio_garcia_garcia/commutativemonoids/
