import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="py4py",
    version="0.9.0",
    author="Python Collaboration",
    author_email="s.mangham@soton.ac.uk",
    description="Analysing and plotting module for the radiative transfer code called Python.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/agnwinds/py4py",
    packages=setuptools.find_packages(),
    scripts=[
        "scripts/tssproduce.py"
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    requires=[
        "numpy",
        "astropy",
        "matplotlib",
        "sqlalchemy",
    ],
    python_requires='>=3.6',
)
