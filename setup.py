from setuptools import setup

setup(
    name="PerturbLoader",
    version="0.0.1",
    py_modules=["Perturb_Loader", ],
    install_requires=["pandas", "anndata", "gdown"],
)
