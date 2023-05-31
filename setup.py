from setuptools import find_packages, setup

setup(
    name="superearth",
    version="1.1",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    package_data={"superearth": ["Data/*.txt"]},
    include_package_data=True,
    python_requires=">=3.6",
    install_requires=[
        "matplotlib",
        "numpy",
        "pandas",
        "scipy"
    ],
    extras_require={
        "optional": [
            "plotly"
        ]
    },
)