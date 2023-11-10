from setuptools import setup

setup(
    name="lrbound",
    version="0.0.1",
    author="Arthur Braida",
    package_dir={"": "src"},
    install_requires=["networkx", "numpy", "qutip", "matplotlib"],
)
