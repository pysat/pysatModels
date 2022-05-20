Contributing
============

Bug reports, feature suggestions and other contributions are greatly
appreciated!  pysatModels is a community-driven project and welcomes both
feedback and contributions.

Short version
=============

* Submit bug reports and feature requests at
  [GitHub Issues](https://github.com/pysat/pysatModels/issues)
* Make pull requests to the ``develop`` branch

Bug reports
===========

When [reporting a bug](https://github.com/pysat/pysatModels/issues) please
include:

* Your operating system name and version
* Any details about your local setup that might be helpful in troubleshooting
* Detailed steps to reproduce the bug

Feature requests and feedback
=============================

The best way to send feedback is to file an issue at
[GitHub Issues](https://github.com/pysat/pysatModels/issues).

If you are proposing a feature:

* Explain in detail how it would work.
* Keep the scope as narrow as possible, to make it easier to implement.
* Remember that this is a volunteer-driven project, and that code contributions
  are welcome :)

Development
===========

To set up `pysatModels` for local development:

- Fork pysat on [GitHub](https://github.com/pysat/pysatModels/fork).

- Clone your fork locally::


       git clone git@github.com:your_name_here/pysatModels.git

- Create a branch for local development::


       git checkout -b name-of-your-bugfix-or-feature

Now you can make your changes locally. Tests for new instruments are
performed automatically.  Tests for custom functions should be added to the
appropriately named file in ``pysatModels/tests``.  For example, the
averaging routines in avg.py are tested in
``pysatModels/tests/test_avg.py``.  If no test file exists, then you should
create one.  This testing uses pytest, which will run tests on any python
file in the test directory that starts with ``test_``.

- When you're done making changes, run all the checks from the
  ``pysatModels/tests`` directory to ensure that nothing is broken on your
  local system.  You may need to install
  [pytest](https://docs.pytest.org/en/latest/) and
  [pytest-flake8](https://pypi.org/project/pytest-flake8/) first. ::


       python -m pytest -vs --flake8

- Update or add documentation (in ``docs``), if relevant.  If you have added
  a new routine, you will need to add an example in the ``docs/examples``
  folder.

- Commit your changes and push your branch to GitHub.  Our commit statements
  follow the basic rules in the
  [Numpy/SciPy workflow](https://docs.scipy.org/doc/numpy-1.15.1/dev/gitwash/development_workflow.html)::


       git add .
       git commit -m "TYPE: Brief description of your changes"
       git push origin name-of-your-bugfix-or-feature

- Submit a pull request through the GitHub website. Pull requests should be
  made to the ``develop`` branch.

Pull Request Guidelines
-----------------------

If you need some code review or feedback while you're developing the code, just
make a pull request. Pull requests should be made to the ``develop`` branch.

For merging, you should:

1. Include an example for use
2. Add a note to ``CHANGELOG.md`` about the changes
3. Ensure that all checks passed (current checks include GitHub Actions
   and Coveralls).

If you don't have all the necessary Python versions available locally or
have trouble building all the testing environments, you can rely on
GitHub to run the tests for each change you add in the pull request.
Because testing here will delay tests by other developers, please ensure
that the code passes all tests on your local system first.
