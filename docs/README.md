# CASTOR-ETC Docs

This is a general read-me for how to get acquainted and set-up with the API documentation. 

## Installation
Required packages for running the documentation website locally is included in the overall project `pyproject.toml`. Please install the `docs` environment using your Python package manager.

> NOTE: it is recommended to do this within a virtual environment so you don't install unnecessary pieces of the 

```Bash
$ pip install -e ".[docs]"
```

## Running locally
After installing required packages, go into the `docs` folder where the `Makefile` exists. You can then create a local version of the wiki with the Makefile.

If you're on Windows:
```Bash
$ ./make.bat html
```

For Unix based systems (macOS / Linux)
```Bash
$ make html
```

This will generate a local copy of the HTML files under the `docs/build/html` folder, which you can then view in a browser. For example, you can run:

```Bash
$ firefox docs/build/html/index.html
```

## Auto-building

If you'd like to have it autobuild on any changes, you can try to set-up [this package](https://github.com/sphinx-doc/sphinx-autobuild)

Just note that the auto-API tooling I use to generate the API documentation *might* not be the easiest to integrate with it.

### Extensions used
> TODO: finish this section!

Full details are found in the [Sphinx configuration Python file](./source/conf.py) found in the `./source/` folder:

- [MySTNB](https://myst-nb.readthedocs.io/en/latest/) : compiles Jupyter notebooks and supports Markdown parsing
