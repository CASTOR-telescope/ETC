GitHub Actions CI/CD
===================

The project uses GitHub Actions for continuous integration and deployment with `uv <https://docs.astral.sh/uv/>`_ for fast dependency management.

Workflows
---------

**Tests** (``test.yml``)
  Runs on every push and pull request. Tests the package across Python 3.10-3.12 on Linux, Windows, and macOS. Generates HTML coverage reports uploaded as artifacts.

**Release** (``release.yml``)
  Runs only on beta, RC, and stable version tags. Runs tests, builds packages with Hatch, and creates GitHub releases.

Creating Releases
-----------------

Create and push a version tag to trigger workflows:

.. code-block:: bash

   # Alpha release (CI testing only, no GitHub release)
   git tag -a v1.0.0-alpha1 -m "Alpha 1" && git push origin v1.0.0-alpha1

   # Beta release (external testing, creates pre-release on GitHub)
   git tag -a v1.0.0-beta1 -m "Beta 1" && git push origin v1.0.0-beta1

   # Release candidate (final testing, creates pre-release on GitHub)
   git tag -a v1.0.0-rc1 -m "RC 1" && git push origin v1.0.0-rc1

   # Stable release (creates full release on GitHub)
   git tag -a v1.0.0 -m "Release 1.0.0" && git push origin v1.0.0

**Release Workflow:**

- **Alpha** → Tests only, no GitHub release (internal CI validation)
- **Beta/RC** → Tests, builds, creates pre-release on GitHub (external testing)
- **Stable** → Tests, builds, creates full release on GitHub

Dependencies
------------

Dependabot automatically checks for dependency updates weekly and creates pull requests for Python packages and GitHub Actions.
