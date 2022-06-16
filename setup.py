from setuptools import setup, find_packages

# Get the version from versioneer
__version__ = "0.0.0"

setup(name="cproofutils",
      version=__version__,
      description="cproof processing utilities",
      author="Jody Klymak",
      author_email="jklymak@gmail.com",
      url="https://github.com/c-proof/cproofutils.git",
      packages=find_packages(exclude=['tests']),
      python_requires='>=3.6',
      install_requires=[
          "pyglider"
      ],
      zip_safe=True
      )
