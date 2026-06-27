import pathlib
import os

cwd = pathlib.Path(__file__).parent

# def pytest_runtest_setup(item):
#     # called for running each test in 'a' directory
#     #os.chdir(cwd)
#     pass


def pytest_configure(config):
    """
    Allows plugins and conftest files to perform initial configuration.
    This hook is called for every plugin and initial conftest
    file after command line options have been parsed.
    """
    os.chdir(cwd)
    print(f"cwd set to {cwd} in {__file__}")


# def pytest_sessionstart(session):
#     """
#     Called after the Session object has been created and
#     before performing collection and entering the run test loop.
#     """


# def pytest_sessionfinish(session, exitstatus):
#     """
#     Called after whole test run finished, right before
#     returning the exit status to the system.
#     """


# def pytest_unconfigure(config):
#     """
#     called before test process is exited.
#     """
