Local Environment Setup
=======================

Introduction
----------------
This document provides instructions for setting up a local development environment for the FORECASTOR ETC tool. It covers the installation of necessary dependencies, configuration of the development environment, and guidelines for contributing to the project.

To start, make sure you have the following prerequisites installed on your system:
- Python 3.9 or higher
- Git
- Python virtual environment tool (e.g., venv or conda)
- A text editor or IDE of your choice (e.g., VSCode, PyCharm)

.. note::
    This document assumes some working knowledge with terminals and Git. While it might seem incredibly frustrating to have to deal with the jargon-filled world of software development, I promise that it is not as hard as it seems.
    `CodeRefinery <https://coderefinery.org/>`_ in general is a great resource to learn scientific computing for researchers. Here are two specific tutorials that will get you up to speed on what "Git" is and how to utilize it for collaborative work.
    - `Introduction to Version Control <https://coderefinery.github.io/git-intro/>`_
    - `Collaborative coding with version control <https://coderefinery.github.io/git-collaborative/>`_

Cloning the Repository
---------------------------
.. note::
    If you are unfamiliar with Git and GitHub, consider going through the `Introduction to Version Control <https://coderefinery.github.io/git-intro/>`_ tutorial from CodeRefinery linked above.

First, you need to clone the FORECASTOR ETC repository from GitHub. Navigate to your working folder and then run the following command in your terminal:

.. code-block:: bash
    git clone https://github.com/CASTOR-telescope/ETC.git

Or, if you prefer GUI-based options, you can use GitHub Desktop or other Git clients to clone the repository.

What cloning the repository does is that it creates a local "instance" of a "Git repository " on your computer. This means that you can make changes to the code, commit those changes, and push them back to the remote repository on GitHub when you're ready.

Side note on Git submodules
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The FORECASTOR ETC repository uses Git submodules to manage example notebooks. This is done to "isolate" the downstream uses from the package itself.
When you clone the repository, you also need to initialize and update the submodules to ensure that all necessary code is available. You can do this by running the following commands in the root directory of the cloned repository:
.. code-block:: bash
    git submodule update --init --recursive

You should then see contents inside the `examples/` folder. This is meant to be a link to a separate repository that contains example notebooks demonstrating the usage of the FORECASTOR ETC. To actually edit these, you will need to clone that repository separately.

Hatchling, our back-end build system
-----------------------------------------------

The FORECASTOR ETC uses `Hatchling <https://hatch.pypa.io/latest/>`_ as its build system and package manager. Hatchling simplifies the process of managing dependencies, building, and distributing Python packages.

In layman's term, this tool will:
- Help you create different "work environments" and switch between them: think of it as moving between notebooks for the same task
- Set the standard and provide a common interface for "building" the package: think of this as creating a version-ed "recipe book" that people can use to replicate the version of the package on their Python

.. note:: 
    If you're curious about the work that went behind setting it up for FORECASTOR ETC, you can check out the official `Python Packaging User Guide <https://packaging.python.org/en/latest/tutorials/packaging-projects/>`_.

To start, install Hatchling on your computer using pip:
.. code-block:: bash
    pip install hatch

Once Hatchling is installed, restart your terminal and verify that your shell can find and run by using the following command:
.. code-block:: bash
    hatch --version
    # Output: 1.12.0

As mentioned above, Hatchling includes features to do "environment management". This is the recommended method of setting up your own environment locally (since I don't want to introduce 3 different tools in this tutorial), but you can feel free to use conda or venv or other tools that you're more familiar with.
A "Python environment" is an isolated runtime that includes a specific Python interpreter, the project's installed packages, and related settings. It prevents dependency and version conflicts between projects and can be created with tools like venv, conda, or Hatch's environment manager.

To create an "environment" using Hatch, enter the ETC root directory and run the following command:
.. code-block:: bash
    hatch env create

Hatch by default creates a virtual environment and also installs all the dependencies that you need. You can "activate" the environment using:
.. code-block:: bash
    hatch shell

This will "open" the environment you just created. Your terminal should now resemble something in the block below. What this essentially does, is point this specific terminal's commands through the virtual environment. Your Python will utilize the environment's Python, which is separate from other Python installations that you might have on the computer.
.. code-block:: bash
    (castor-etc) $

`castor_etc` will be installed in the environment, and you can verify this by running:
.. code-block:: bash
    pip show castor_etc

Manually creating a virtual environment
---------------------------------------

If you prefer to manually build the virtual environment, you can do so using Python's built-in venv module or conda. Here are the steps for both methods:

Using venv
1. Navigate to the root directory of the cloned FORECASTOR ETC repository.
2. Create a virtual environment by running:
.. code-block:: bash
    python -m venv castor-etc

3. Activate the virtual environment:
.. code-block:: bash

    # On Windows
    .\castor-etc\Scripts\activate

    # On macOS/Linux
    source castor-etc/bin/activate

4. Install required dependencies by 
.. code-block:: bash

    pip install -e .

Using conda
1. Navigate to the root directory of the cloned FORECASTOR ETC repository.
2. Create a conda environment by running:
.. code-block:: bash
    conda create --name castor-etc python=3.x
3. Activate the conda environment:
.. code-block:: bash
    conda activate castor-etc
4. Install required dependencies by
.. code-block:: bash
    pip install -e .


