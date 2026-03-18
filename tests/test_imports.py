import pkgutil
import castor_etc
import pytest

def test_package_imports():
    """ Recursively check that all submodules can be imported. """
    package = castor_etc
    for loader, module_name, is_pkg in pkgutil.walk_packages(package.__path__, package.__name__ + "."):
        try:
            __import__(module_name)
        except Exception as e:
            pytest.fail(f"Failed to import {module_name}: {e}")