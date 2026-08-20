Writing & Running Automated Tests in Python
===========================================

The CASTOR ETC uses pytest for testing. This short guide is here to help you get set-up with the test environment locally and understand how to debug it.

Folder Structure
----------------

Python automated tests are located in the root of the ETC repository, under the **tests** folder.

The folder contains a few package level tests along with a folder for unit tests (**unit**). For the most part, most developers should only be concerned with the higher level tests. Unit testing is meant for checking individual components.


Adding New Tests
----------------

To add a new test, create a file under **tests** with a filename that starts with ``test__``.

Inside the file, you can either use the `unittest <https://docs.python.org/3/library/unittest.html>`_ or `PyTest <https://docs.pytest.org/en/stable/>`_ subpackages, which both provide detailed guides on how to write your test files.

You can also reference the `point source photometry test file <https://github.com/CASTOR-telescope/ETC/blob/master/tests/test_point_source_photometry.py>`_ for how existing tets are set-up.

Running Tests
-------------

Assuming you have `Hatch <https://hatch.pypa.io/latest/>`_ installed, you can run the following command in your terminal:

.. code-block:: bash

    $ hatch test

If you don't have Hatch installed and would prefer to do things the old-school way, run these commands in your terminal:

.. code-block:: bash
    
    $ python -m pip install -e ".[test]"
    $ pytest tests/