===========
MFLike-camspec
===========

An ``MFLike`` implementation of the Planck ``camspec-npipe`` high-ell likelihood for `cobaya <https://github.com/CobayaSampler/cobaya>`_.

Running/testing the code
------------------------

You can test the ``MFLike-CamSpec`` likelihood by including the ``CamSpecMFLike`` in your ``likelihoods`` block. All parameters used by ``CamSpec`` have the same name in ``MFLike-CamSpec``.

You need to unzip the ``data/like_NPIPE_12.6_unified_cov.bin.zip`` file to include the ``data/like_NPIPE_12.6_unified_cov.bin`` file. After that you can test the code by running

.. code:: shell

    $ python3 -m cobaya run mflike-camspec.yaml

Included files
--------------

The following files and folders are included:

    ``camspecmflike/`` contains the main ``CamSpecMFLike`` code base.

    ``camspecmflike/fgspectra/`` contains a custom ``fgspectra`` implementation that is needed for the correct ``camspec`` foreground modeling. It contains the Planck foreground models in the ``fgspectra/data/`` folder. Foreground model implementations are included that replicate the Planck foregrounds based on the original templates.

    ``data/`` is the default path to place the ``camspec-NPIPE`` data files (you can also download them from the `CamSpec 2020 Public Release <https://people.ast.cam.ac.uk/~stg20/camspec/index.html>`_).

    ``ini_yamls/`` contains several ``.yaml`` files that replicate the chains ran for this paper, while ``yamls/`` contains the set up defaults for these chains.

Two example/test ``.yaml`` files are included:

    ``mflike-camspec.yaml`` runs a basic mcmc chain over the data.
    
    ``mflike-camspec-lowLE.yaml`` couples the ``MFLike-CamSpec`` code to the ``plik-lowl`` and ``plik-lowE`` codes and performs a combined chain over them. It should be modified to point to the local ``.clik`` files (see the `cobaya documentation <https://cobaya.readthedocs.io/en/latest/likelihood_planck.html>`_ for more information on that).
    
    The remaining ``params_*.yaml`` files define the parameters and priors as originally presented in the `cobaya NPIPE highl CamSpec <https://github.com/CobayaSampler/cobaya/tree/master/cobaya/likelihoods/planck_NPIPE_highl_CamSpec>`_ implementation.
    
The sample chains should be able to be ran as

.. code:: shell

    $ python3 -m cobaya run mflike-camspec.yaml

It should be noted that the code should be run using ``python 3``, it is incompatible with ``python 2``. For the ``mflike-camspec-lowLE.yaml`` example, you should direct it to your local ``.clik`` files as mentioned above.
