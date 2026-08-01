import pathlib
import os

cwd = pathlib.Path(__file__).parent


def pytest_configure(config):
    """
    Allows plugins and conftest files to perform initial configuration.
    This hook is called for every plugin and initial conftest
    file after command line options have been parsed.
    """
    os.chdir(cwd)
    print(f"cwd set to {cwd} in {__file__}")
