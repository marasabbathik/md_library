from setuptools import setup, find_packages

from setuptools import setup, find_packages

setup(
    name="md_library",
    version="0.1.0",
    description="Physics data analysis library for XRD, RSM, ellipsometry and visualization",
    author="Marek Duchan",
    author_email="mduchan@mail.muni.cz",
    url="https://github.com/yourusername/md_library",  # Optional
    packages=find_packages(),
    python_requires=">=3.8",
    install_requires=[
        "numpy>=1.20.0",
        "matplotlib>=3.3.0",
        # Add other dependencies your library needs
    ],
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Physics",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
    ],
)