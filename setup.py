from setuptools import setup, find_packages
import codecs
import os

here = os.path.abspath(os.path.dirname(__file__))

version_file = open(os.path.join(here, 'VERSION'))
version = version_file.read().strip()

# Get the long description from the relevant file
with codecs.open(os.path.join(here, 'DESCRIPTION.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='pynucl',

    # Versions should comply with PEP440. For single-sourced versioning, see
    # http://packaging.python.org/en/latest/tutorial.html#version
    version=version,

    description='A library for comprehensive analysis of nucleosome structures',
    long_description=long_description,

    # The project URL.
    url='https://github.com/intbio/pynucl',

    # Author details
    author='Alexey K. Shaytan',
    author_email='alex@intbio.org',

    # Choose your license
    license='GPL 3.0',

    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        # https://pypi.python.org/pypi?%3Aaction=list_classifiers
        # Project maturity. 
        'Development Status :: 3 - Alpha',

        # Intended audience
        'Intended Audience :: Science/Research',

        # Topic
        'Topic :: Scientific/Engineering :: Bio-Informatics',

        # License should match "license" above
        'License :: GPL 3.0',

        # Python versions supported
        'Programming Language :: Python :: 3',
    ],

    # What does your project relate to?
    keywords='science sequences bioinformatics',

    # You can just specify the packages manually here if your project is
    # simple. Or you can use find_packages().
    packages=find_packages(exclude=['examples', 'docs', 'tests*']),
    # Run-time package dependencies. These will be installed by pip when your
    # project is installed.
    install_requires=[
        'biopython'
    ],

    # Data files included in your packages. If using Python 2.6 or less, 
    # then these have to be included in MANIFEST.in as well.
    include_package_data=True,
    # To provide executable scripts, use entry points in preference to the
    # "scripts" keyword. Entry points provide cross-platform support and allow
    # pip to create the appropriate form of executable for the target platform.
    python_requires='==3.*.*',

    # Default to False unless you specifically intend the package to be
    # installed as an .egg
    zip_safe=False,
)
