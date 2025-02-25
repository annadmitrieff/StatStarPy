from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="statstar",
    version="0.1.0",
    author="Anna Dmitrieff",
    author_email="annadmitrieff@uga.edu",
    description="Python translation of the Fortran StatStar stellar modeling code, from 'An Introduction to Modern Astrophysics.'",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/annadmitrieff/statstar",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
    install_requires=[
        "numpy",
        "scipy",
        "matplotlib",
    ],
)