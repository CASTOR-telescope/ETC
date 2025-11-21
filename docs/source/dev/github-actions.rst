GitHub Actions CI/CD
===================

The project uses GitHub Actions for continuous integration and deployment.

Workflows
---------

**Tests** (``test.yml``)
  Runs on every push and pull request. Tests the package across Python 3.9-3.12 on Linux, Windows, and macOS. Generates coverage reports and posts them as PR comments.

**Build Check** (``build.yml``)
  Runs only on version tags. Verifies the package builds correctly before release.

**Release** (``release.yml``)
  Runs only on version tags. Builds multi-platform wheels, creates GitHub releases, and optionally publishes to PyPI.

Creating Releases
-----------------

Create and push a version tag to trigger a release:

.. code-block:: bash

   # Alpha release (for testing builds)
   git tag -a v1.0.0-alpha -m "Alpha 1" && git push origin v1.0.0-alpha1

   # Beta release (feature-complete testing)
   git tag -a v1.0.0-beta1 -m "Beta 1" && git push origin v1.0.0-beta1

   # Stable release
   git tag -a v1.0.0 -m "Release 1.0.0" && git push origin v1.0.0

Alpha and beta releases are automatically marked as pre-releases on GitHub.

Dependencies
------------

Dependabot automatically checks for dependency updates weekly and creates pull requests for Python packages and GitHub Actions.
